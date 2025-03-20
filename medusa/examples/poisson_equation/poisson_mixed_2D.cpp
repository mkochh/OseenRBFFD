#include <medusa/Medusa_fwd.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

/// Basic medusa example, we are solving 2D Poisson's equation on unit square with
/// mixed boundary conditions in order to calculate
/// quarter of the solution to the Dirichlet Boundary conditions (because of symmetry).
/// Neumann boundary conditions on the upper and right side of the box
/// and Dirichlet boundary conditions on the bottom and left side of the box
/// http://e6.ijs.si/medusa/wiki/index.php/Poisson%27s_equation

using namespace mm;  // NOLINT

int main() {
    // Create the domain and discretize it
    BoxShape<Vec2d> box(0.0, 0.5);  // Square with a = 0.5
    double dx = 0.005;
    DomainDiscretization<Vec2d> domain = box.discretizeWithStep(dx);

    // Find support for the nodes
    int N = domain.size();
    domain.findSupport(FindClosest(9));  // the support for each node is the closest 9 nodes

    // Construct the approximation engine, in this case a weighted least squares using Gaussian RBF,
    // no weight, scale to farthest and LLT solver.
    WLS<Gaussians<Vec2d>, NoWeight<Vec2d>, ScaleToFarthest,
            Eigen::LLT<Eigen::MatrixXd>> wls({9, 30});

    // Laplacian for inside the domain and first derivative for the bc's
    auto storage = domain.computeShapes<sh::lap|sh::d1>(wls);

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    Eigen::VectorXd rhs(N); rhs.setZero();
    M.reserve(storage.supportSizes());

    auto op = storage.implicitOperators(M, rhs);  // construct implicit operators over our storage

    int BOTTOM = -3;
    int TOP = -4;
    int LEFT = -1;
    int RIGHT = -2;

    for (int i : domain.interior()) {
        double x = domain.pos(i, 0);
        double y = domain.pos(i, 1);
        op.lap(i) = -2*PI*PI*std::sin(PI*x)*std::sin(PI*y);  // set the case for nodes in the domain
    }
    for (int i : (domain.types() == LEFT) + (domain.types() == BOTTOM)) {
        op.value(i) = 0.0;  // Dirichlet boundary conditions on the left and bottom edge of the box
    }
    for (int i : (domain.types() == TOP) + (domain.types() == RIGHT)) {
        // Neumann boundary conditions on upper and right edge of the box
        op.neumann(i, domain.normal(i)) = 0.0;
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);

    std::ofstream out_file("poisson_mixed_2D_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "solution = " << u << ";" << std::endl;
    out_file.close();

    return 0;
}
