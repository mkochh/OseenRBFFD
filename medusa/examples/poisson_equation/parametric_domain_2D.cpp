#include <medusa/Medusa.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

/// Basic medusa example, we are solving 2D Poisson's equation on a parametric domain
/// discretized with variable density with Dirichlet boundary conditions.
/// http://e6.ijs.si/medusa/wiki/index.php/Parametric_domains

using namespace mm;  // NOLINT

int main() {
    // Define paramteric curve and its Jacobian matrix.
    auto example_r = [](Vec<double, 1> t) {
        double r = pow(abs(cos(1.5 * t(0))), sin(3 * t(0)));
        return Vec2d(r * cos(t(0)), r * sin(t(0)));
    };

    auto der_example_r = [](Vec<double, 1> t) {
        double r = pow(abs(cos(1.5 * t(0))), sin(3 * t(0)));
        double der_r = (-1.5 * pow(abs(cos(1.5 * t(0))),
                sin(3 * t(0))) * sin(3 * t(0)) * sin(1.5 * t(0)) +
                        3 * pow(abs(cos(1.5 * t(0))),
                                sin(3 * t(0))) * cos(3 * t(0)) * cos(1.5 * t(0))
                                * log(abs(cos(1.5 * t(0))))) / cos(1.5 * t(0));

        Eigen::Matrix<double, 2, 1> jm;
        jm.col(0) << der_r * cos(t(0)) - r * sin(t(0)), der_r * sin(t(0)) + r * cos(t(0));

        return jm;
    };

    // Define parametric curve's domain.
    BoxShape<Vec1d> param_bs(Vec<double, 1>{0.0}, Vec<double, 1>{2 * PI});

    // Fill domain.
    UnknownShape<Vec2d> shape;
    DomainDiscretization<Vec2d> domain(shape);

    auto gradient_h = [](Vec2d p){
        double h_0 = 0.005;
        double h_m = 0.03 - h_0;

        return  (0.5 * h_m * (p(0) + p(1) + 3.0) + h_0) / 5.0;
    };

    GeneralSurfaceFill<Vec2d, Vec1d> gsf;
    domain.fill(gsf, param_bs, example_r, der_example_r, gradient_h);

    GeneralFill<Vec2d> gf;
    domain.fill(gf, gradient_h);

    // Find support for the nodes.
    int N = domain.size();
    domain.findSupport(FindClosest(9));  // the support for each node is the closest 9 nodes

    // Construct the approximation engine.
    int m = 2;  // basis order
    Monomials<Vec2d> mon(m);
    RBFFD<Polyharmonic<double, 3>, Vec2d> approx({}, mon);

    // Compute the shapes (we only need the Laplacian).
    auto storage = domain.computeShapes<sh::lap>(approx);

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    Eigen::VectorXd rhs(N); rhs.setZero();
    M.reserve(storage.supportSizes());

    // Construct implicit operators over our storage.
    auto op = storage.implicitOperators(M, rhs);

    for (int i : domain.interior()) {
        // set the case for nodes in the domain
        op.lap(i) = 0.5;
    }
    for (int i : domain.boundary()) {
        // enforce the boundary conditions
        op.value(i) = 0.0;
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);

    // Write the solution into file.
    std::ofstream out_file("parametric_domain_2D_data.m");
    out_file << "positions = " << domain.positions() << ";" << std::endl;
    out_file << "solution = " << u << ";" << std::endl;
    out_file.close();
    return 0;
}
