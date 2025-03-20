#include <medusa/Medusa_fwd.hpp>
#include <medusa/bits/domains/BasicRelax.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

/// Basic medusa example, we are solving 2D Poisson's equation on unit circle
/// with mixed boundary conditions
/// http://e6.ijs.si/medusa/wiki/index.php/Poisson%27s_equation

using namespace mm;  // NOLINT

int main() {
    // construct BallShape and discretize the boundary
    BallShape<Vec2d> c(0.0, 1.0);
    double dx = 0.05;
    DomainDiscretization<Vec2d> domain = c.discretizeBoundaryWithStep(dx);

    // fill using GeneralFill, when filling the domain like this
    // it needs to be relaxed afterwards
    GeneralFill<Vec2d> fill;
    domain.fill(fill, dx);

    // relax the domain
    BasicRelax relax;
    relax.iterations(20).initialHeat(0.8).numNeighbours(3)
            .projectionType(BasicRelax::DO_NOT_PROJECT);
    relax(domain, dx);

    int N = domain.size();
    domain.findSupport(FindClosest(9));

    WLS<Gaussians<Vec2d>, NoWeight<Vec2d>, ScaleToFarthest,
            Eigen::LLT<Eigen::MatrixXd>> wls({9, 30});

    auto storage = domain.computeShapes<sh::lap|sh::d1>(wls);

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    Eigen::VectorXd rhs(N); rhs.setZero();

    auto op = storage.implicitOperators(M, rhs);
    M.reserve(storage.supportSizes());

    for (int i : domain.interior()) {
        op.lap(i) = 1.0;
    }
    for (int i : domain.boundary()) {
        double x = domain.pos(i, 0);
        if (x > 0) {
            op.value(i) = 0.0;  // dirichlet
        } else {
            op.neumann(i, domain.normal(i)) = 0.0;  // neumann
        }
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);

    std::ofstream out_file("poisson_mixed_2D_circular_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "solution = " << u << ";" << std::endl;
    out_file.close();

    return 0;
}

