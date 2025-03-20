#include <medusa/Medusa_fwd.hpp>
#include <medusa/bits/domains/BasicRelax.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

/// Basic medusa example, we are solving 2D Poisson's equation
/// on unit square with mixed boundary conditions
/// Using RBFFD, Gaussian RBF augmented with monomials
/// http://e6.ijs.si/medusa/wiki/index.php/Poisson%27s_equation

using namespace mm;  // NOLINT

int main() {
    // Create the domain and discretize it
    BoxShape<Vec2d> b(0.0, 1.0);
    double dx = 0.05;

    DomainDiscretization <Vec2d> domain = b.discretizeBoundaryWithStep(dx);

    // fill the domain
    GeneralFill <Vec2d> fill;
    fill.seed(999);
    domain.fill(fill, dx);

    // relax the domain
    BasicRelax relax;
    relax.iterations(20).initialHeat(0.8).finalHeat(0.0).numNeighbours(3).projectionType(
            BasicRelax::DO_NOT_PROJECT);
    relax(domain, dx);

    int N = domain.size();
    FindClosest find_support(12);
    domain.findSupport(find_support);

    Gaussian<double> g(30);  // construct Gaussians
    // Augmented Gaussian RBF, with monomials up to order 2
    RBFFD<Gaussian<double>, Vec2d, ScaleToClosest> appr(g, Monomials<Vec2d>(2));

    auto storage = domain.computeShapes<sh::lap | sh::d1>(appr);
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    Eigen::VectorXd rhs(N); rhs.setZero();
    Range<int> reserve = storage.supportSizes();
    auto op = storage.implicitOperators(M, rhs);
    M.reserve(storage.supportSizes());

    int BOTTOM = -3;
    int TOP = -4;
    int LEFT = -1;
    int RIGHT = -2;

    for (int i : domain.interior()) {
        double x = domain.pos(i, 0);
        double y = domain.pos(i, 1);
        op.lap(i) = -std::cos(x + y) * std::exp(-x - y);
    }
    for (int i : (domain.types() == LEFT)) {
        double y = domain.pos(i, 1);
        op.neumann(i, domain.normal(i)) = -std::cos(PI * y);
    }
    for (int i : (domain.types() == RIGHT)) {
        double y = domain.pos(i, 1);
        op.neumann(i, domain.normal(i)) = -std::sin(PI * y);
    }
    for (int i : (domain.types() == BOTTOM)) {
        double x = domain.pos(i, 0);
        op.value(i) = -std::exp(x);
    }
    for (int i : (domain.types() == TOP)) {
        op.value(i) = 0.0;
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);

    // Write the solution into file
    std::ofstream out_file("poisson_mixed_2D_RBFFD_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "solution = " << u << ";" << std::endl;
    out_file.close();

    return 0;
}
