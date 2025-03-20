#include <medusa/Medusa_fwd.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

/// Basic medusa example, we are solving 1D Poisson's equation on unit line with Dirichlet
/// boundary conditions.
/// http://e6.ijs.si/medusa/wiki/index.php/Poisson%27s_equation

using namespace mm;  // NOLINT

int main() {
    // Create the domain and discretize it
    BoxShape <Vec1d> box(0.0, 1.0);
    double dx = 0.01;
    DomainDiscretization <Vec1d> domain = box.discretizeWithStep(dx);

    // Find support for the nodes
    int N = domain.size();
    domain.findSupport(FindClosest(4));  // the support for each node is the closest 4 nodes

    // Construct the approximation engine, in this case a weighted least squares using monomials
    // as basis functions, Gaussian weight, and scale to farthest
    int m = 2;  // basis order
    Monomials <Vec1d> basis(m);  // construct monomial basis
    WLS <Monomials<Vec1d>, GaussianWeight<Vec1d>,
    ScaleToFarthest> wls(basis);  // basis of monomials up to order 2

    // compute the shapes (only Laplacian in this case) using our WLS
    auto storage = domain.computeShapes<sh::lap>(wls);

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    Eigen::VectorXd rhs(N); rhs.setZero();
    M.reserve(storage.supportSizes());

    // construct implicit operators over our storage
    auto op = storage.implicitOperators(M, rhs);

    for (int i : domain.interior()) {
        double x = domain.pos(i, 0);
        // set the case for nodes in the domain
        op.lap(i) = -PI * PI * std::sin(PI * x);
    }
    for (int i : domain.boundary()) {
        // enforce the boundary conditions
        op.value(i) = 0.0;
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.compute(M);
    // solve the system
    ScalarFieldd u = solver.solve(rhs);

    // Write the solution into file
    std::ofstream out_file("poisson_dirichlet_1D_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "solution = " << u << ";" << std::endl;
    out_file.close();

    return 0;
}
