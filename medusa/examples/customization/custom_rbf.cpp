#include <medusa/Medusa.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>

/// Custom RBF example. Solves Poisson equation.
/// http://e6.ijs.si/medusa/wiki/index.php/Customization

using namespace mm;  // NOLINT

/// f(r) = r^3 + r^5.
class MyRBF {
  public:
    typedef double scalar_t;
    double operator()(double r2) const {
        double r = std::sqrt(r2);
        return ipow<3>(r) + ipow<5>(r);
    }

    template <int dim>
    double operator()(double r2, Lap<dim>) const {
        double r = std::sqrt(r2);
        return 3*(dim+1)*r + 5*(dim+3)*ipow<3>(r);
    }
};

int main() {
    // Create the domain and discretize it
    BallShape<Vec2d> b(0.0, 1.0);
    double dx = 0.05;

    DomainDiscretization <Vec2d> domain = b.discretizeBoundaryWithStep(dx);

    // fill the domain
    GeneralFill <Vec2d> fill;
    fill.seed(999);
    domain.fill(fill, dx);

    int N = domain.size();
    FindClosest find_support(12);
    domain.findSupport(find_support);

    RBFFD<MyRBF, Vec2d, ScaleToClosest> appr({}, Monomials<Vec2d>(2));

    auto storage = domain.computeShapes<sh::lap>(appr);
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    Eigen::VectorXd rhs(N); rhs.setZero();
    Range<int> reserve = storage.supportSizes();
    auto op = storage.implicitOperators(M, rhs);
    M.reserve(storage.supportSizes());

    for (int i : domain.interior()) {
        double x = domain.pos(i, 0);
        double y = domain.pos(i, 1);
        op.lap(i) = -2*PI*PI*std::sin(x)*std::sin(y);
    }
    for (int i : domain.boundary()) {
        double x = domain.pos(i, 0);
        double y = domain.pos(i, 1);
        op.value(i) = std::sin(x) * std::sin(y);
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);

    // Write the solution into file
    std::ofstream out_file("custom_rbf_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "solution = " << u << ";" << std::endl;
    out_file.close();

    return 0;
}
