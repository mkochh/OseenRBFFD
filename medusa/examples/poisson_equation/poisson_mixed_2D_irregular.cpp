#include <medusa/Medusa_fwd.hpp>
#include <medusa/bits/domains/BasicRelax.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

/// Basic medusa example, we are solving 2D Poisson's equation on polygon domain with
/// 3 edges having dirichlet bc's and 2 having Neumann bc's.
/// http://e6.ijs.si/medusa/wiki/index.php/Poisson%27s_equation

using namespace mm;  // NOLINT

int main() {
    // Create a polygon shape, with 5 points
    PolygonShape <Vec2d> poly({{0.0,  0.0},
                               {2.8,  1.8},
                               {2.3,  4.1},
                               {0.7,  2.5},
                               {-2.3, 1.7}});

    double dx = 0.05;
    // discretize the boundary
    DomainDiscretization <Vec2d> domain = poly.discretizeBoundaryWithStep(dx);

    // fill the inside using GeneralFill
    GeneralFill <Vec2d> fill;
    fill.proximityTolerance(0.9).seed(420);
    domain.fill(fill, dx);

    // relax the domain
    BasicRelax relax;
    relax.initialHeat(0.5).finalHeat(0.1).numNeighbours(3)
            .projectionType(BasicRelax::DO_NOT_PROJECT);
    domain.relax(relax, dx);

    // Find support for the nodes
    int N = domain.size();
    domain.findSupport(FindClosest(9));

    // Construct the approximation engine, in this case a weighted least squares using Gaussian RBF,
    // no weight, scale to farthest and LLT solver.
    WLS<Gaussians<Vec2d>, NoWeight<Vec2d>, ScaleToFarthest, Eigen::LLT<Eigen::MatrixXd>>
    wls({9, 30.0});

    auto storage = domain.computeShapes<sh::lap | sh::d1>(wls);  // compute Lap and d1 shapes

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    Eigen::VectorXd rhs(N); rhs.setZero();

    auto op = storage.implicitOperators(M, rhs);  // construct implicit operators over our storage
    M.reserve(storage.supportSizes());
    for (int i : domain.interior()) {
        op.lap(i) = 1;
    }
    for (int i : domain.boundary()) {
        op.value(i) = 0.0;
    }
    // edges will have negative integer types
    for (int i : (domain.types() == -5) + (domain.types() == -2)) {
        op.neumann(i, domain.normal(i)) = 0.0;
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);

    std::ofstream out_file("poisson_mixed_2D_irregular_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "solution = " << u << ";" << std::endl;
    out_file.close();
    return 0;
}
