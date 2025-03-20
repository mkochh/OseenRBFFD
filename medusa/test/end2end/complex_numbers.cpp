#include <medusa/Medusa_fwd.hpp>
#include <complex>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

#include "gtest/gtest.h"

constexpr std::complex<double> operator""_i(long double d) {
    return std::complex<double>{0.0, static_cast<double>(d)};
}

namespace mm {

TEST(End2end, ComplexNumbers) {
    BoxShape<Vec2d> box(0.0, 1.0);
    DomainDiscretization<Vec2d> domain = box.discretizeWithStep(0.05);

    int N = domain.size();
    domain.findSupport(FindClosest(9));

    WLS<Gaussians<Vec2d>> wls({9, 30.0}, {});
    auto storage = domain.computeShapes(wls);

    // Implicit scalar
    Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> M(N, N);
    M.reserve(storage.supportSizes());
    Eigen::VectorXcd rhs(N); rhs.setZero();

    ImplicitOperators<decltype(storage), decltype(M), decltype(rhs)> op(storage, M, rhs);
    double x, y;
    for (int i : domain.interior()) {
        x = domain.pos(i, 0); y = domain.pos(i, 1);
        1.0_i * op.lap(i) = 2*PI*PI*std::sin(PI * x)*std::sin(PI * y);
    }
    for (int i : domain.boundary()) { op.value(i) = 0.0; }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<std::complex<double>>> solver;
    solver.compute(M);
    Eigen::VectorXcd u = solver.solve(rhs);

    EXPECT_NEAR(u.real().norm(), 0, 1e-10);
    for (int i : domain.interior()) {
        x = domain.pos(i, 0); y = domain.pos(i, 1);
        double e = std::abs(u[i] - 1.0_i*std::sin(PI * x)*std::sin(PI * y));
        EXPECT_NEAR(e, 0, 1.5e-2);
    }

    // Implicit vector
    M.setZero();
    M.reserve(storage.supportSizes());
    rhs.setZero();
    ImplicitOperators<decltype(storage), decltype(M), decltype(rhs)> opv(storage, M, rhs);
    for (int i : domain.interior()) {
        x = domain.pos(i, 0); y = domain.pos(i, 1);
        1.0_i * opv.lap(i) = 2*PI*PI*std::sin(PI * x)*std::sin(PI * y);
    }
    for (int i : domain.boundary()) { opv.value(i) = 0.0; }

    solver.compute(M);
    Eigen::VectorXcd uv = solver.solve(rhs);
    ASSERT_EQ(uv, u);

    auto vf = VectorField<std::complex<double>, 1>::fromLinear(uv);

    auto eop = storage.explicitOperators();
    auto evop = storage.explicitVectorOperators();
    for (int i : domain.interior()) {
        x = domain.pos(i, 0); y = domain.pos(i, 1);
        // Explicit scalar
        auto egrad = eop.grad(u, i);
        auto agrad = Eigen::Vector2cd(1.0_i*PI*std::cos(PI*x)*std::sin(PI*y),
                                      1.0_i*PI*std::sin(PI*x)*std::cos(PI*y));
        double e = (egrad - agrad).norm();
        EXPECT_NEAR(e, 0, 1e-2);

        // Explicit vector
        auto m = evop.grad(vf, i).transpose();
        EXPECT_EQ(m, egrad);
    }
}

}   // namespace mm
