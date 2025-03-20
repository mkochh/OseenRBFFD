#include <medusa/Medusa_fwd.hpp>
#include <Eigen/SparseCore>
#include <Eigen/Geometry>
#include <Eigen/IterativeLinearSolvers>

/// Basic medusa example, we are solving 3D Poisson's equation on rotated and translated unit cube
/// with Dirichlet boundary conditions.
/// http://e6.ijs.si/medusa/wiki/index.php/Poisson%27s_equation

using namespace mm;  // NOLINT

int main() {
    // Create the domain and discretize it
    BoxShape <Vec3d> box(0.0, 1.0);
    double dx = 0.05;
    DomainDiscretization <Vec3d> domain = box.discretizeWithStep(dx);

    Eigen::AngleAxisd Q(PI/3, Vec3d(1, 1, 1).normalized());
    domain.rotate(Q.toRotationMatrix()).translate({-2.5, 0.4, 1.2});

    // Find support for the nodes
    int N = domain.size();
    domain.findSupport(FindClosest(9));  // the support for each node is the closest 9 nodes

    // Construct the approximation engine, in this case a weighted least squares using Gaussians
    // as basis functions, no weight, and scale to farthest
    WLS <Gaussians<Vec3d>, NoWeight<Vec3d>, ScaleToFarthest,
    Eigen::LLT<Eigen::MatrixXd>> wls({9, 30.0});

    auto storage = domain.computeShapes<sh::lap>(wls);  // compute the shapes using our WLS

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    Eigen::VectorXd rhs(N); rhs.setZero();

    auto op = storage.implicitOperators(M, rhs);  // construct implicit operators over our storage
    M.reserve(storage.supportSizes());

    for (int i : domain.interior()) {
        double x = domain.pos(i, 0);
        double y = domain.pos(i, 1);
        double z = domain.pos(i, 2);
        // set the case for nodes in the domain
        op.lap(i) = -3 * PI * PI * std::sin(PI * x) * std::sin(PI * y) * std::sin(PI * z);
    }
    for (int i : domain.boundary()) {
        // enforce the boundary conditions
        op.value(i) = 0.0;
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);  // solve the system

    // Write the solution into file
    std::ofstream out_file("poisson_dirichlet_3D_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "solution = " << u << ";" << std::endl;
    out_file.close();
    return 0;
}

