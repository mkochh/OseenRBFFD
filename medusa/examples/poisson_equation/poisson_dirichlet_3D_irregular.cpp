#include <medusa/Medusa_fwd.hpp>
#include <medusa/bits/domains/BasicRelax.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

/// Solution of the Poisson equation on an irregular domain in 3D, with mixed bc's.
/// This example will combine all the knowledge from the previous ones along with
/// presenting XML and HDF classes
/// that are key for I/O in bigger programs.

using namespace mm;  // NOLINT

int main() {
    XML conf("poisson_dirichlet_3D_irregular_params.xml");

    // Obtain parameters
    double dx = conf.get<double>("num.dx");
    int n = conf.get<int>("approximations.n");
    int m = conf.get<int>("approximations.m");
    double sigma = conf.get<double>("approximations.sigma");

    BoxShape<Vec3d> box(0.0, 1.0);
    BallShape<Vec3d> s1({0.8, 0.8, 0.9}, 0.5);  // sphere to subtract from the box

    DomainDiscretization <Vec3d> domain = box.discretizeBoundaryWithStep(dx);
    auto d = s1.discretizeBoundaryWithStep(dx);
    domain -= d;  // subtract the domains

    // fill the interior using GeneralFill
    GeneralFill <Vec3d> fill;
    domain.fill(fill, dx);
    BasicRelax relax;

    // relax the domain
    relax.iterations(20).initialHeat(0.8).numNeighbours(3).projectionType(
            BasicRelax::DO_NOT_PROJECT);
    relax(domain, dx);

    int N = domain.size();
    domain.findSupport(FindClosest(n));

    WLS<Gaussians<Vec3d>, NoWeight<Vec3d>, ScaleToFarthest,
    Eigen::LLT<Eigen::MatrixXd>> wls({m, sigma});

    auto storage = domain.computeShapes<sh::lap>(wls);

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    Eigen::VectorXd rhs(N); rhs.setZero();

    auto op = storage.implicitOperators(M, rhs);
    M.reserve(storage.supportSizes());

    for (int i : domain.interior()) {
        double x = domain.pos(i, 0);
        double y = domain.pos(i, 1);
        double z = domain.pos(i, 2);
        op.lap(i) = -3 * PI * PI * std::sin(PI * x) * std::sin(PI * y) * std::sin(PI * z);
    }
    for (int i : domain.boundary()) {
        op.value(i) = 0.0;
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);

    HDF hdf("poisson_dirichlet_3D_irregular.h5", HDF::DESTROY);
    hdf.writeDoubleArray("solution", u);
    hdf.writeDouble2DArray("positions", domain.positions());
    hdf.close();

    return 0;
}
