#include <medusa/Medusa_fwd.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

/// The study of fretting fatigue, for which the authors believe no closed form solution
/// is known.  A small thin rectangular specimen of width $W$, length $L$ and
/// thickness $t$ is stretched in one axis with axial traction $\sigma_{ax}$ and
/// compressed in another with force $F$ by two oscillating cylindrical pads of
/// radius $R$, inducing tangential force $Q$. The tractions induced by the pads are predicted
/// using an extension of Hertzian contact theory, which splits the contact area into
/// the stick and slip zones, depending on the coefficient of friction $\mu$ and
/// combined Young's modulus $E^*$, given by $\frac{1}{E^*} = \frac{1 - \nu_1^2}{E_1} +
/// \frac{1 - \nu_2^2}{E_2}$, where $E_i$ and $\nu_i$ represent the
/// Young's moduli and Poisson's ratios of the specimen and the pad, respectively.
///
/// Example taken from: https://onlinelibrary.wiley.com/doi/full/10.1002/nme.6067

using namespace mm;  // NOLINT

int main() {
    // Load params.
    double E = 72.1e9;
    double nu = 0.33;
    double F = 543;
    double Q = 155;
    double COF = 0.3;
    double thickness = 0.004;
    double radius = 0.01;
    double sigma_axial = 1e8;

    // Parameter logic.
    double mu = E / 2. / (1 + nu);
    double lam = E * nu / (1 - 2 * nu) / (1 + nu);
    double Estar = E / (2 * (1 - nu * nu));
    double a = 2 * std::sqrt(std::abs(F * radius / (thickness * PI * Estar)));
    double p0 = std::sqrt(std::abs(F * Estar / (thickness * PI * radius)));
    double c = a * std::sqrt(1 - Q / COF / std::abs(F));
    double e = a * sigma_axial / 4 / COF / p0;

    // Domain definition.
    const double width = 0.04;
    const double height = 0.01;
    const double dx = 0.00005;
    BoxShape<Vec2d> shape({-width / 2, -height / 2}, {width / 2, 0});

    auto h = [=](const Vec2d &) { return dx; };
    DomainDiscretization<Vec2d> domain = shape.discretizeBoundaryWithDensity(h);

    GeneralFill<Vec2d> fill_randomized;
    int seed = 100;
    fill_randomized.seed(seed);
    fill_randomized(domain, h);

    int N = domain.size();
    prn(N);
    domain.findSupport(FindClosest(12));  // The support for each node is the closest 12 nodes.

    // Construct the approximation engine, in this case a RBF-FD with PHS and monomials.
    RBFFD<Polyharmonic<double>, Vec2d, ScaleToClosest> engine(3, 2);

    // Compute the shapes using the defined approximation engine
    auto storage = domain.computeShapes(engine);
    std::cerr << "Shapes computed, solving system." << std::endl;

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(2 * N, 2 * N);
    Eigen::VectorXd rhs(2 * N);
    rhs.setZero();
    M.reserve(storage.supportSizesVec());

    // Construct implicit operators over our storage.
    auto op = storage.implicitVectorOperators(M, rhs);

    // Set the governing equations and the boundary conditions.
    for (int i : domain.interior()) {
        (lam + mu) * op.graddiv(i) + mu * op.lap(i) = 0.0;
    }
    for (int i : domain.boundary()) {
        int id = domain.type(i);
        double x = domain.pos(i, 0);
        double trac = 0;

        switch (id) {
            case -4:
                // TOP.
                op.traction(i, lam, mu, {0, 1}) = 1.0;
                if (c <= std::abs(x + e) && std::abs(x) <= a) {
                    trac = -COF * p0 * std::sqrt(1 - (x / a) * (x / a));
                } else if (std::abs(x + e) < c) {
                    trac = -COF * p0 *
                           (std::sqrt(1 - (x / a) * (x / a)) -
                            c / a * std::sqrt(1 - (x + e) * (x + e) / c / c));
                }
                rhs(i) = trac;
                rhs(i + N) = (std::abs(x) < a) ? -p0 * std::sqrt(1 - x * x / a / a) : 0;
                break;
            case -2:
                // RIGHT.
                op.traction(i, lam, mu, {1, 0}) = 1.0;
                rhs(i) = sigma_axial;
                rhs(i + N) = 0;
                break;
            case -3:
                // BOTTOM.
                op.eq(0).c(0).der1(i, 1) = 0.0;
                op.eq(1).c(1).value(i) = 0;
                rhs(i) = 0;
                M.coeffRef(i + N, i + N) = 1;
                rhs(i + N) = 0;
                break;
            default:
                // LEFT.
                op.value(i) = {0.0, 0.0};
                rhs(i) = rhs(i + N) = 0;
        }
    }

    Eigen::VectorXd norms = (M.cwiseAbs() * Eigen::VectorXd::Ones(M.cols())).cwiseInverse();
    M = norms.asDiagonal() * M;
    rhs = rhs.cwiseProduct(norms);

    // Prepare solver.
    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.preconditioner().setFillfactor(40);
    solver.preconditioner().setDroptol(1e-5);
    solver.setMaxIterations(100);
    solver.setTolerance(1e-15);
    solver.compute(M);

    // Solve.
    Eigen::VectorXd sol = solver.solve(rhs);
    prn(solver.iterations());
    prn(solver.error());

    // Postprocess.
    VectorField2d u = VectorField2d::fromLinear(sol);
    VectorField3d stress(N);
    auto eop = storage.explicitVectorOperators();
    for (int i = 0; i < N; ++i) {
        auto grad = eop.grad(u, i);
        stress(i, 0) = (2 * mu + lam) * grad(0, 0) + lam * grad(1, 1);
        stress(i, 1) = lam * grad(0, 0) + (2 * mu + lam) * grad(1, 1);
        stress(i, 2) = mu * (grad(0, 1) + grad(1, 0));
    }

    // Write the solution into a file.
    HDF hdf("results_fatigue.h5", HDF::DESTROY);
    hdf.writeDomain("domain", domain);
    hdf.writeEigen("displacements", u);
    hdf.writeEigen("stress", stress);
    hdf.closeFile();

    return 0;
}
