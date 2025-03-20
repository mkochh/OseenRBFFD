#include <medusa/Medusa_fwd.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/LU>

/// Basic medusa example, we are solving 2D Poisson's equation on unit square
/// with Dirichlet boundary conditions.
/// http://e6.ijs.si/medusa/wiki/index.php/Poisson%27s_equation

using namespace mm;  // NOLINT

int main() {
    // Create the domain and discretize it
    BallShape<Vec2d> box(0.0, 1.0);
    double dx = 0.05;
    DomainDiscretization <Vec2d> domain = box.discretizeWithStep(dx);

    for (int level = 1; level <= 3; ++level) {
        domain.findSupport(FindClosest(7));
        HalfLinksRefine refine;
        refine.region(domain.positions().filter([=](const Vec2d& v) {
            return v.norm() < 1.0 / (level+1);
        }));
        refine(domain);
    }

    // Find support for the nodes
    int N = domain.size();
    domain.findSupport(FindClosest(13));  // the support for each node is the closest 9 nodes

    // Construct the approximation engine, in this case a weighted least squares using monomials as
    // basis functions, no weight, and scale to farthest
    int m = 2;  // basis order
    RBFFD<Polyharmonic<double, 3>, Vec2d, ScaleToClosest> approx({}, m);

    // compute the shapes (we only need the Laplacian) using our WLS
    auto storage = domain.computeShapes<sh::lap>(approx);

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    Eigen::VectorXd rhs(N); rhs.setZero();

    // construct implicit operators over our storage
    auto op = storage.implicitOperators(M, rhs);
    M.reserve(storage.supportSizes());
    double s = 0.1;
    for (int i : domain.interior()) {
        double r = domain.pos(i).norm();
        -op.lap(i) = std::exp(-r*r/s/s);
    }
    for (int i : domain.boundary()) {
        op.value(i) = 0.0;
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);

    // Write the solution into file
    std::ofstream out_file("poisson_dirichlet_2D_refined_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "solution = " << u << ";" << std::endl;
    out_file.close();
    return 0;
}
