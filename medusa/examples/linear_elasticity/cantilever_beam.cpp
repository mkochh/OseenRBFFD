#include <medusa/Medusa_fwd.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>


/// Classical 2D Cantilever beam problem.
/// http://e6.ijs.si/medusa/wiki/index.php/Linear_elasticity#Cantilever_beam

using namespace mm;  // NOLINT

int main() {
    // Physical parameters
    const double E = 72.1e9;
    const double nu = 0.33;
    const double P = 1000;
    const double D = 5;
    const double L = 30;

    // Derived parameters
    const double I = D * D * D / 12;
    double lam = E * nu / (1 - 2 * nu) / (1 + nu);
    const double mu = E / 2 / (1 + nu);
    lam = 2 * mu * lam / (2 * mu + lam);  // plane stress

    // Edge labels
    int LEFT = -1;
    int RIGHT = -2;
    int BOTTOM = -3;
    int TOP = -4;

    // Domain definition
    Vec2d low(0, -D / 2), high(L, D / 2);
    BoxShape<Vec2d> box(low, high);

    double step = 0.25;
    DomainDiscretization<Vec2d> domain = box.discretizeWithStep(step);
    int N = domain.size();
    prn(N);
    domain.findSupport(FindClosest(9));  // the support for each node is the closest 9 nodes

    // Construct the approximation engine, in this case a weighted least squares using monomials as
    // basis functions, no weight, and scale to farthest.
    int m = 2;  // basis order
    WLS<Monomials<Vec2d>, NoWeight<Vec2d>, ScaleToFarthest>
    wls(Monomials<Vec2d>::tensorBasis(m));

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
    for (int i : domain.types() == RIGHT) {
        double y = domain.pos(i, 1);
        op.value(i) = {(P*y*(3*D*D*(1+nu) - 4*(2+nu)*y*y)) / (24.*E*I), -(L*nu*P*y*y) / (2.*E*I)};
    }
    for (int i : domain.types() == LEFT) {
        double y = domain.pos(i, 1);
        op.traction(i, lam, mu, {-1, 0}) = {0, -P*(D*D - 4*y*y) / (8.*I)};
    }
    for (int i : domain.types() == TOP) {
        op.traction(i, lam, mu, {0, 1}) = 0.0;
    }
    for (int i : domain.types() == BOTTOM) {
        op.traction(i, lam, mu, {0, -1}) = 0.0;
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
    std::ofstream out_file("cantilever_beam_data.m");
    out_file << "L = " << L << ";\n";
    out_file << "D = " << D << ";\n";
    out_file << "E = " << E << ";\n";
    out_file << "nu = " << nu << ";\n";
    out_file << "P = " << P << ";\n";
    out_file << "positions = " << domain.positions() << ";\n";
    out_file << "displacement = " << u << ";" << ";\n";
    out_file << "stress = " << stress << ";" << std::endl;
    out_file.close();
    return 0;
}
