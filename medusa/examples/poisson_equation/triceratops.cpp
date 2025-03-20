#include <medusa/Medusa_fwd.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

/// Medusa example of loading a pre-prepared domain that can be imported and solved with medusa.
/// As long as a hdf5 file is available that contains the positions, types, normals and a
/// boundary map between normals and types this can be done.
/// http://e6.ijs.si/medusa/wiki/index.php/Poisson%27s_equation

using namespace mm;  // NOLINT

int main() {
    // Import domain from h5 file
    HDF hdf("triceratops_domain.h5");
    DomainDiscretization<Vec3d> domain = DomainDiscretization<Vec3d>::load(hdf, "/domain");
    hdf.close();

    // Find support for the nodes
    int N = domain.size();
    domain.findSupport(FindClosest(15));  // the support for each node is the closest 9 nodes

    // Construct the approximation engine, in this case a weighted least squares using Gaussians
    // as basis functions, no weight, and scale to farthest
    WLS<Gaussians<Vec3d>, NoWeight<Vec3d>, ScaleToFarthest> wls({15, 50.0});

    auto storage = domain.computeShapes<sh::lap>(wls);  // compute the shapes using our WLS

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    Eigen::VectorXd rhs(N); rhs.setZero();

    auto op = storage.implicitOperators(M, rhs);  // construct implicit operators over our storage
    M.reserve(storage.supportSizes());
    for (int i : domain.interior()) {
        op.lap(i) = -1;
    }
    for (int i : domain.boundary()) {
        op.value(i) = 0.0;
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);  // solve the system

    // Write the solution into file
    HDF out("triceratops.h5", HDF::DESTROY);
    out.writeDoubleArray("solution", u);
    out.close();

    return 0;
}
