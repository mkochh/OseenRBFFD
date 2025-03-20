#include <medusa/Medusa.hpp>
#include <medusa/bits/domains/PolyhedronShape.hpp>  // This must be included separately.
#include <Eigen/Sparse>

// Example using PolyhedronShape with realistic 3D objects.
// http://e6.ijs.si/medusa/wiki/index.php/Realistic_3D_models
// Needs CGAL (https://www.cgal.org/).

using namespace mm;  // NOLINT

int main() {
    std::string filename = "../../test/testdata/bunny.off";
    PolytopeShape<Vec3d> bunny = PolytopeShape<Vec3d>::fromOFF(filename);

    double dx = 2.5;
    DomainDiscretization<Vec3d> domain = bunny.discretizeBoundaryWithStep(dx);
    GeneralFill<Vec3d> fill;
    fill.seed(0);
    fill(domain, dx);

    int N = domain.size();
    prn(N);
    domain.findSupport(FindClosest(25));
    RBFFD<Polyharmonic<double, 3>, Vec3d, ScaleToClosest> approx({}, 2);
    auto storage = domain.template computeShapes<sh::lap|sh::d1>(approx);

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    Eigen::VectorXd rhs(N); rhs.setZero();

    auto op = storage.implicitOperators(M, rhs);
    M.reserve(storage.supportSizes());

    for (int i : domain.interior()) {
        0.1*op.grad(i, -1) + 2.0*op.lap(i) = -1.0;
    }
    for (int i : domain.boundary()) {
        op.value(i) = 0.0;
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.preconditioner().setDroptol(1e-4);
    solver.preconditioner().setFillfactor(20);
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);

    HDF hdf("bunny_poisson.h5", HDF::DESTROY);
    hdf.writeDouble2DArray("pos", domain.positions());
    hdf.writeEigen("sol", u);

    return 0;
}
