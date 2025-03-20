#include <medusa/Medusa_fwd.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>


/// Classical 2D Point contact problem.
/// http://e6.ijs.si/medusa/wiki/index.php/Linear_elasticity#Point_contact_2D

using namespace mm;  // NOLINT

int main() {
    // Physical parameters
    const double E = 72.1e9;
    const double nu = 0.33;
    const double P = 1000;

    // Derived parameters
    double lam = E * nu / (1 - 2 * nu) / (1 + nu);
    const double mu = E / 2 / (1 + nu);

    // Closed form solution
    std::function<Vec2d(Vec2d)> analytical = [=](const Vec2d& p) {
        double x = p[0], y = p[1], r2 = x*x + y*y, factor = -P / (4.0*PI*mu);
        double ux = 2*x*y/r2 + 2*mu/(lam + mu)*std::atan2(y, x);
        double uy = (y*y-x*x)/r2 - (lam + 2*mu)/(lam + mu)*std::log(r2);
        return Vec2d(factor*ux, factor*uy);
    };

    // Domain definition
    Vec2d low(-1, -1), high(1, -0.1);
    BoxShape<Vec2d> box(low, high);

    double step = 0.025;
    DomainDiscretization<Vec2d> domain = box.discretizeBoundaryWithStep(step);
    GeneralFill<Vec2d> fill; fill.seed(0);
    fill(domain, step);
    int N = domain.size();
    prn(N);
    domain.findSupport(FindClosest(9));  // the support for each node is the closest 9 nodes

    // Construct the approximation engine, in this case a weighted least squares using monomials as
    // basis functions, no weight, and scale to farthest.
    WLS<Gaussians<Vec2d>, GaussianWeight<Vec2d>, ScaleToFarthest> wls({9, 100.0}, 1.0);

    // Compute the shapes using the defined approximation engine.
    auto shapes = domain.computeShapes(wls);
    std::cerr << "Shapes computed, solving system." << std::endl;

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(2*N, 2*N);
    Eigen::VectorXd rhs(2*N); rhs.setZero();
    M.reserve(shapes.supportSizesVec());

    // Construct implicit operators over our storage.
    auto op = shapes.implicitVectorOperators(M, rhs);

    // Set the governing equations and the boundary conditions.
    for (int i : domain.interior()) {
        (lam+mu)*op.graddiv(i) + mu*op.lap(i) = 0.0;
    }
    for (int i : domain.boundary()) {
        op.value(i) = analytical(domain.pos(i));
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.preconditioner().setFillfactor(20);
    solver.preconditioner().setDroptol(1e-4);
    solver.compute(M);
    Eigen::VectorXd solution = solver.solve(rhs);
    prn(solver.iterations());
    prn(solver.error());
    VectorField2d u = VectorField2d::fromLinear(solution);

    VectorField3d stress(N);
    auto eop = shapes.explicitVectorOperators();
    for (int i = 0; i < N; ++i) {
        auto grad = eop.grad(u, i);
        stress(i, 0) = (2*mu + lam)*grad(0, 0) + lam*grad(1, 1);
        stress(i, 1) = lam*grad(0, 0) + (2*mu+lam)*grad(1, 1);
        stress(i, 2) = mu*(grad(0, 1)+grad(1, 0));
    }

    // Write the solution into a file.
    std::ofstream out_file("point_contact2d_data.m");
    out_file << "E = " << E << ";\n";
    out_file << "nu = " << nu << ";\n";
    out_file << "P = " << P << ";\n";
    out_file << "positions = " << domain.positions() << ";\n";
    out_file << "displacement = " << u << ";" << ";\n";
    out_file << "stress = " << stress << ";" << std::endl;
    out_file.close();
    return 0;
}
