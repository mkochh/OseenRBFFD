#include <medusa/Medusa_fwd.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <complex>

// Electromagnetic scattering on cylinder example.
// http://e6.ijs.si/medusa/wiki/index.php/Electromagnetic_scattering

using namespace mm;  // NOLINT

constexpr std::complex<double> operator"" _i(long double d) {
    return std::complex<double>{0.0, static_cast<double>(d)};
}

int main() {
    double dx = 0.01;
    double m0 = 1.256 * 1e-6;  // magnetic permeability
    double e0 = 8.85 * 1e-12;
    double f = 3 * 1e8;  // linear frequency
    double ww = 2 * PI * f;  // angular frequency
    double er = 0.5;  // specific permittivity
    double t0 = PI/4;  // incident angle
    // A matrix entries
    double mxx = 2;
    double myy = 2;
    double mxy = 0.75;
    double myx = 0.75;
    double ga = (mxx * myy - mxy * myx);
    mxx = mxx / ga;
    myy = myy / ga;
    mxy = mxy / ga;
    myx = myx / ga;

    double ka = PI/2;  // size factor
    double tol = 1e-5;
    double k0 = ww * sqrt(m0 * e0);  // free space wave vector
    double r1 = ka / k0;  // inner radius
    double r2 = 3*r1;  // outer radius

    int n = 9;

    // domain
    BallShape<Vec2d> innerc(0, r1);
    BallShape<Vec2d> outerc(0, r2);
    auto annulus = outerc - innerc;

    // discretize annulus boundary
    DomainDiscretization<Vec2d> outer = annulus.discretizeBoundaryWithStep(dx);

    // Filter inner boundary nodes
    DomainDiscretization<Vec2d> inner(innerc);
    auto inner_bnd = outer.positions().filter(
            [=](const Vec2d& p) { return std::abs(p.norm() - r1) < tol; });
    auto outer_bnd = outer.positions().filter(
            [=](const Vec2d& p) { return std::abs(p.norm() - r2) < tol; });

    // Index mapping from outer to inner
    Range<int> outer_to_inner(outer.size());
    for (int i = 0; i < static_cast<int>(inner_bnd.size()); ++i) {
        outer_to_inner[inner_bnd[i]] = i;
        inner.addBoundaryNode(outer.pos(inner_bnd[i]), -2, -1 * outer.normal(inner_bnd[i]));
    }

    GeneralFill<Vec2d> fill_engine;
    inner.fill(fill_engine, dx);  // the annulus
    outer.fill(fill_engine, dx);  // the inner circle

    // Find support
    inner.findSupport(FindClosest(n));
    outer.findSupport(FindClosest(n));

    // Create WLS
    WLS<Gaussians<Vec2d>, NoWeight<Vec2d>, ScaleToFarthest> wls({n, 30});

    // Create complex sparse matrix
    int N_inner = inner.size();
    int N_outer = outer.size();

    // Create matrix and rhs
    Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> M(N_inner + N_outer,
                                                                 N_inner + N_outer);
    Eigen::VectorXcd rhs(N_inner + N_outer);

    // shape storage for both domains
    auto storage_inner = inner.computeShapes<sh::lap | sh::d1 | sh::d2>(wls);
    auto storage_outer = outer.computeShapes<sh::lap | sh::d1 | sh::d2>(wls);

    // operators for both domains
    auto op_inner = storage_inner.implicitOperators(M, rhs);
    auto op_outer = storage_outer.implicitOperators(M, rhs);

    op_inner.setRowOffset(N_outer);
    op_inner.setColOffset(N_outer);

    // m*2 in each row, for
    auto store = storage_inner.supportSizes() + storage_outer.supportSizes();
    Range<int> reserve = store;
    for (int& i : reserve) i *= 2;

    M.reserve(reserve);
    for (int i : inner.interior()) {
        // Wave eq on the inner circle
        // (mxx*v_xx + myy*v_yy + mxy*v_xy + myx*v_yx) + k0^2*er*v = 0
        mxx * op_inner.der2(i, 0, 0) + myy * op_inner.der2(i, 1, 1) +
        (mxy + myx) * op_inner.der2(i, 0, 1)
        + k0 * k0 * er * op_inner.value(i) = 0;
    }
    for (int i : inner.boundary()) {
        // dirichlet, u_inside - u_scattered = u_incident<- incoming wave from t0 angle
        double x = inner.pos(i, 0);
        double y = inner.pos(i, 1);
        double theta = atan2(y, x);  // get phase

        // in the upper side of the matrix
        std::complex<double> inc = std::exp(1.0_i * k0 * (x * std::cos(t0) + y * std::sin(t0)));
        op_inner.value(i) + (-1) * op_outer.value(inner_bnd[i], N_outer + i) = inc;

        // Neumann boundary condition on the inside circle
        auto norm = inner.normal(i);
        double n0 = norm[0];
        double n1 = norm[1];

        // Anisotropic normal derivative
        // (n0*mxx+n1*myx)*v_x + (n0*mxy+n1*myy)*v_y + du/dn = d(incident)/dn
        op_outer.neumann(inner_bnd[i], outer.normal(inner_bnd[i]))
        + (n0 * mxx + n1 * myx) * op_inner.der1(i, 0, -N_outer + inner_bnd[i])
        + (n0 * mxy + n1 * myy) * op_inner.der1(i, 1, -N_outer + inner_bnd[i])
                = 1.0_i * k0 * std::cos(theta - t0) *
                  std::exp(1.0_i * k0 * (x * std::cos(t0) + y * std::sin(t0)));
    }
    for (int i : outer.interior()) {  // wave equation for the outer region
        op_outer.lap(i) + k0 * k0 * op_outer.value(i) = 0;
    }
    for (int i : outer_bnd) {  // Sommerfeld boundary condition
        op_outer.neumann(i, outer.normal(i)) + (k0 * 1.0_i + 1.0/(2*r2)) * op_outer.value(i) = 0.0;
    }

    HDF out("anisotropic_cylinder.h5", HDF::DESTROY);
    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<std::complex<double>>> solver;
    solver.compute(M);

    Eigen::VectorXcd sol = solver.solve(rhs);
    std::cout << "Iterations: " << solver.iterations() << std::endl;
    std::cout << "Error: " << solver.error() << std::endl;

    ScalarField<double> rsol = sol.real();
    ScalarField<double> csol = sol.imag();

    // Configuration data
    out.writeDoubleAttribute("dx", dx);
    out.writeDoubleAttribute("a", r1);
    Eigen::SparseMatrix<double, Eigen::RowMajor> R(N_inner + N_outer, N_inner + N_outer);

    // Data for Domain Plot
    out.writeDoubleArray("rsol", rsol);
    out.writeDoubleArray("csol", csol);
    out.writeDomain("inner", inner);
    out.writeDomain("outer", outer);
    out.close();

    return 0;
}
