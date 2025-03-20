#include <medusa/Medusa_fwd.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>


/// Basic medusa example, we are solving 2D Poisson's equation on unit square
/// with Dirichlet boundary conditions.
/// http://e6.ijs.si/medusa/wiki/index.php/Poisson%27s_equation

using namespace mm;  // NOLINT

int main() {
    // Create the domain and discretize it
    BoxShape <Vec2d> box(0.0, 1.0);
    double dx = 0.01;
    DomainDiscretization <Vec2d> domain = box.discretizeWithStep(dx);

    // Find support for the nodes
    int N = domain.size();
    domain.findSupport(FindClosest(9));  // the support for each node is the closest 9 nodes

    // Construct the approximation engine, in this case a weighted least squares using monomials as
    // basis functions, no weight, and scale to farthest
    int m = 2;  // basis order
    WLS <Monomials<Vec2d>, NoWeight<Vec2d>,
    ScaleToFarthest> wls(Monomials<Vec2d>::tensorBasis(m));  // tensor basis of monomials
    // up to order 2,
    // {1,x,x2,y,yx,yx2,y2,y2x,y2x2}

    // compute the shapes (we only need the Laplacian) using our WLS
    auto storage = domain.computeShapes<sh::lap>(wls);

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    Eigen::VectorXd rhs(N); rhs.setZero();
    M.reserve(storage.supportSizes());

    // construct implicit operators over our storage
    auto op = storage.implicitOperators(M, rhs);

    for (int i : domain.interior()) {
        double x = domain.pos(i, 0);
        double y = domain.pos(i, 1);
        // set the case for nodes in the domain
        op.lap(i) = -2 * PI * PI * std::sin(PI * x) * std::sin(PI * y);
    }
    for (int i : domain.boundary()) {
        // enforce the boundary conditions
        op.value(i) = 0.0;
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);

    // Write the solution into file
    std::ofstream out_file("poisson_dirichlet_2D_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "solution = " << u << ";" << std::endl;
    out_file.close();
    return 0;
}
