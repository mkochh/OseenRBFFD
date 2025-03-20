#include <medusa/Medusa_fwd.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

/// Basic medusa example, that solves the Poisson equation with only Neumann boundary conditions.
/// When the compatibility condition is met (integral of the source function inside the domain is
/// equal to the net flow expressed by the boundary integral of the normal derivative), we need
/// another constraint in order to obtain a unique solution
/// This is done by regularization (Average of the field = 0). Details:
/// http://e6.ijs.si/medusa/wiki/index.php/Poisson%27s_equation

using namespace mm;  // NOLINT

int main() {
    BoxShape <Vec2d> b(0.0, 1.0);
    double dx = 0.05;
    DomainDiscretization <Vec2d> domain = b.discretizeWithStep(dx);

    // Find support for the nodes
    int N = domain.size();
    domain.findSupport(FindClosest(12));  // the support for each node is the closest 9 nodes

    // Construct the approximation engine, in this case a weighted least squares using Gaussian RBF,
    // no weight, scale to farthest
    Monomials <Vec2d> basis(2);
    WLS <Monomials<Vec2d>, GaussianWeight<Vec2d>, ScaleToFarthest> wls(basis);

    // Laplacian an d1 shape computation
    auto storage = domain.computeShapes<sh::lap | sh::d1>(wls);

    // extra row and column for the lagrange multiplier
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N + 1, N + 1);

    // extra row for the lambda parameter
    Eigen::VectorXd rhs(N+1); rhs.setZero();

    // construct implicit operators over our storage
    auto op = storage.implicitOperators(M, rhs);
    for (int i : domain.interior()) {
        op.lap(i) = 2;
    }

    int BOTTOM = -3;
    int TOP = -4;
    int LEFT = -1;
    int RIGHT = -2;

    for (int i : (domain.types() == LEFT)) {
        op.der1(i, 0) = 0.0;
    }
    for (int i : (domain.types() == RIGHT)) {
        op.der1(i, 0) = 1.0;
    }
    for (int i : (domain.types() == BOTTOM)) {
        op.der1(i, 1) = 0.0;
    }
    for (int i : (domain.types() == TOP)) {
        op.der1(i, 1) = 1.0;
    }
    // regularization
    // set the last row and column of the matrix
    for (int i = 0; i < N; ++i) {
        M.coeffRef(N, i) = 1;
        M.coeffRef(i, N) = 1;
    }
    // set the sum of all values
    rhs[N] = 0.0;

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.compute(M);
    ScalarFieldd solution = solver.solve(rhs);
    ScalarFieldd u = solution.head(N);

    // Write the solution into file
    std::ofstream out_file("poisson_neumann_2D_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "solution = " << u << ";" << std::endl;
    out_file << "lambda = " << solution(N) << ";" << std::endl;
    out_file.close();

    return 0;
}

