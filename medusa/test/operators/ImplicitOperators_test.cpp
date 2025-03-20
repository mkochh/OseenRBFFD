#include <medusa/bits/approximations/JacobiSVDWrapper_fwd.hpp>
#include <medusa/bits/approximations/Monomials.hpp>
#include <medusa/bits/approximations/ScaleFunction.hpp>
#include <medusa/bits/approximations/WeightFunction_fwd.hpp>
#include <medusa/bits/approximations/Gaussian_fwd.hpp>
#include <medusa/bits/approximations/WLS.hpp>
#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/domains/FindClosest.hpp>
#include <medusa/bits/operators/computeShapes.hpp>
#include <medusa/bits/operators/ImplicitOperators.hpp>
#include <medusa/bits/operators/RaggedShapeStorage.hpp>
#include <medusa/bits/operators/UniformShapeStorage.hpp>
#include <medusa/bits/types/Vec_fwd.hpp>
#include <medusa/bits/operators/ShapeStorage.hpp>


#include <Eigen/LU>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

#include <functional>
#include <cmath>

#include "gtest/gtest.h"

namespace mm {

class ImplicitOp : public ::testing::Test {
  public:
    UniformShapeStorage<Vec2d> shapes;
    Eigen::MatrixXd M;
    Eigen::VectorXd rhs;
    ImplicitOperators<decltype(shapes), decltype(M), decltype(rhs)> op;

    ImplicitOp() : M(3, 3), rhs(3) {}

  protected:
    void SetUp() override {
        shapes.resize({3, 3, 3});
        shapes.setSupport(0, {0, 1, 2});
        shapes.setSupport(1, {0, 1, 2});
        shapes.setSupport(2, {0, 1, 2});

        M.setZero();
        rhs.setZero();

        Eigen::Vector3d s0(1.2, 1.3, 1.4), s1(-14, -15, -16), s2(100, 110, 120);
        shapes.setLaplace(0, s0);
        shapes.setLaplace(1, s1);
        shapes.setLaplace(2, s2);

        Eigen::Vector3d d0(1200, 1300, 1400), d1(-14000, -15000, -16000), d2(0.01, 0.011, 0.012);
        shapes.setD1(0, 0, d0);
        shapes.setD1(0, 1, d1);
        shapes.setD1(0, 2, d2);
        shapes.setD1(1, 0, 2*d0);
        shapes.setD1(1, 1, 2*d1);
        shapes.setD1(1, 2, 2*d2);

        Eigen::Vector3d dd0(0.5, 0.6, 0.7), dd1(-1.4, -1.5, -1.6), dd2(-25, -26, -27);
        shapes.setD2(0, 0, 0, dd0);
        shapes.setD2(0, 0, 1, dd1);
        shapes.setD2(0, 0, 2, dd2);
        shapes.setD2(0, 1, 0, 2*dd0);
        shapes.setD2(0, 1, 1, 2*dd1);
        shapes.setD2(0, 1, 2, 2*dd2);
        shapes.setD2(1, 1, 0, 3*dd0);
        shapes.setD2(1, 1, 1, 3*dd1);
        shapes.setD2(1, 1, 2, 3*dd2);

        op = shapes.implicitOperators(M, rhs);
    }
};

TEST_F(ImplicitOp, AssignTwice) {
    Eigen::Matrix3d expectedM = M;
    Eigen::Vector3d expectedRhs = rhs;

    expectedM(0, 0) = 1;
    expectedRhs[0] = 1;
    op.value(0) = 1;
    EXPECT_EQ(0, (expectedM - M).lpNorm<1>());
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());

    expectedM(0, 0) = 2;
    expectedRhs[0] = 2;
    op.value(0) = 1;
    EXPECT_EQ(0, (expectedM - M).lpNorm<1>());
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());
}

TEST_F(ImplicitOp, AssignCombined) {
    double v = 4.2;
    Eigen::Matrix3d expectedM = M;
    Eigen::Vector3d expectedRhs = rhs;
    op.value(1) + v * op.lap(1) + op.value(1) = v;
    expectedM.row(1) += v * shapes.laplace(1);
    expectedM(1, 1) += 2;
    expectedRhs[1] = v;
    EXPECT_EQ(0, (expectedM - M).lpNorm<1>());
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());
}

TEST_F(ImplicitOp, AssignValue) {
    double v = 4.2, f = -3.4;
    Eigen::Matrix3d expectedM = M;
    Eigen::Vector3d expectedRhs = rhs;
    expectedM(1, 1) = f;
    expectedRhs[1] = v;
    f*op.value(1) = v;
    EXPECT_EQ(0, (expectedM - M).lpNorm<1>());
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());
}

TEST_F(ImplicitOp, AssignLap) {
    double v = 4.2;
    Eigen::Matrix3d expectedM = M;
    Eigen::Vector3d expectedRhs = rhs;
    expectedM.row(0) = shapes.laplace(0);
    expectedRhs[0] = v;
    op.lap(0) = v;
    EXPECT_EQ(0, (expectedM - M).lpNorm<1>());
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());
}

TEST_F(ImplicitOp, AssignGrad) {
    double v = 4.2;
    Vec2d direction(3.4, 2.4);
    Eigen::Matrix3d expectedM = M;
    Eigen::Vector3d expectedRhs = rhs;
    expectedM.row(0) = direction[0]*shapes.d1(0, 0) + direction[1]*shapes.d1(1, 0);
    expectedRhs[0] = v;
    op.grad(0, direction) = v;
    EXPECT_EQ(0, (expectedM - M).lpNorm<1>());
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());
}

TEST_F(ImplicitOp, AssignD1) {
    double v = 4.2;
    Eigen::Matrix3d expectedM = M;
    Eigen::Vector3d expectedRhs = rhs;
    expectedM.row(0) = shapes.d1(0, 0);  // dx in node 0
    expectedM.row(1) = shapes.d1(1, 1);  // dy in node 1
    expectedRhs[0] = v;
    expectedRhs[1] = v;
    op.der1(0, 0) = v;
    op.der1(1, 1) = v;
    EXPECT_EQ(0, (expectedM - M).lpNorm<1>());
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());
}

TEST_F(ImplicitOp, AssignD2) {
    double v = 4.2;
    Eigen::Matrix3d expectedM = M;
    Eigen::Vector3d expectedRhs = rhs;
    expectedM.row(0) = shapes.d2(0, 0, 0);  // dxdx in node 0
    expectedM.row(1) = shapes.d2(0, 1, 1);  // dxdy in node 1
    expectedM.row(2) = shapes.d2(1, 1, 2);  // dydy in node 2
    expectedRhs[0] = v;
    expectedRhs[1] = v;
    expectedRhs[2] = v;
    op.der2(0, 0, 0) = v;
    op.der2(1, 0, 1) = v;
    op.der2(2, 1, 1) = v;
    EXPECT_EQ(0, (expectedM - M).lpNorm<1>());
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());
}

TEST_F(ImplicitOp, AssignDeath) {
    EXPECT_DEATH(op.lap(0) + op.lap(1) = 5, "Cannot add together terms for different matrix rows");
}

TEST_F(ImplicitOp, Syntax) {
    constexpr int n = 7;
    Eigen::MatrixXd newM(n, 3); newM.setZero();
    Eigen::VectorXd newRhs(n); newRhs.setZero();

    op.setProblem(newM, newRhs, 0, 0);
    op.lap(0) = 1;
    1.0*op.lap(0, 1) = 1;
    0.0*op.lap(0, 2) + op.lap(0, 2) = 1;
    op.lap(0, 3) + 0.0*op.lap(0, 3) = 1;
    0.25*op.lap(0, 4) + 0.75*op.lap(0, 4) = 1;
    -op.lap(0, 5) + 2.0*op.lap(0, 5) = 1;
    2.0*op.lap(0, 6) + -op.lap(0, 6) = 1;

    for (int i = 1; i < n; ++i) {
        EXPECT_EQ(newM.row(i), newM.row(0));
        EXPECT_EQ(newRhs.row(i), newRhs.row(0));
    }
}

TEST(Operators, ImplicitComplexNumbers) {
    UniformShapeStorage<Vec1d> shapes;
    shapes.resize({2, 2});

    shapes.setSupport(0, {0, 1});
    shapes.setSupport(1, {0, 1});

    Eigen::Vector2d s0(-0.12, 2.42), s1(0.34, 23.3);
    shapes.setLaplace(0, s0);
    shapes.setLaplace(1, s1);

    Eigen::Matrix2cd M; M.setZero();
    Eigen::Vector2cd rhs; rhs.setZero();

    auto op = shapes.implicitOperators(M, rhs);
    std::complex<double> c1(1.2, -3.4), c2(0.7, 8.9), c3(-5.6, 2.3);
    double x = 1.6;
    c1*op.lap(0) + c2*op.value(0) = c3;
    op.lap(1) + c1*op.lap(1) + x*op.lap(1) = c3;

    Eigen::Matrix2cd expectedM = M;
    Eigen::Vector2cd expectedRhs = rhs;
    expectedM.row(0) = c1*s0;
    expectedM(0, 0) += c2;
    expectedRhs[0] = c3;

    expectedM.row(1) = (c1+1.0+x)*s1;
    expectedRhs[1] = c3;

    EXPECT_LT((expectedM - M).lpNorm<1>(), 3e-16);
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());
}

TEST(Operators, Implicit1DLaplaceDense) {
    /// [Implicit laplace 1d example]
    // Solve 1-D boundary value problem
    // u''(x) = x; u(0) = 1, u(1) = -2
    // with solution u(x) = 1/6 (6 - 19 x + x^3)
    BoxShape<Vec1d> box(0.0, 1.0);
    DomainDiscretization<Vec1d> domain = box.discretizeWithStep(0.1);
    WLS<Monomials<Vec1d>> wls(4);
    int N = domain.size();
    domain.findSupport(FindClosest(5));
    auto shapes = domain.computeShapes<sh::lap>(wls, domain.interior());

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(N, N);
    Eigen::VectorXd rhs(N); rhs.setZero();

    auto op = shapes.implicitOperators(M, rhs);

    for (int i : domain.interior()) {
        op.lap(i) = domain.pos(i, 0);
    }
    for (int i : domain.boundary()) {
        op.value(i) = (domain.pos(i, 0) == 0) ? 1 : -2;
    }
    Eigen::VectorXd sol = M.lu().solve(rhs);
    std::function<double(double)> analytic = [](double x) { return 1/6.0 * (6 - 19 * x + x*x*x); };
    for (int i = 0; i < domain.size(); ++i) {
        EXPECT_NEAR(sol[i], analytic(domain.pos(i, 0)), 1e-14);
    }
    /// [Implicit laplace 1d example]
}

TEST(Operators, Implicit2dLaplaceSparse) {
    // Solve boundary value problem
    // laplace u = 0, u on edge = 1
    // with solution u == 1
    BoxShape<Vec2d> box(0.0, 1.0);
    DomainDiscretization<Vec2d> domain = box.discretizeWithStep(0.2);
    WLS<Monomials<Vec2d>> wls(2);
    int N = domain.size();
    domain.findSupport(FindClosest(9));
    auto shapes = domain.computeShapes<sh::lap>(wls, domain.interior());

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    M.reserve(shapes.supportSizes());
    Eigen::VectorXd rhs(N); rhs.setZero();

    auto op = shapes.implicitOperators(M, rhs);

    /// [Implicit syntax]
    for (int i : domain.interior()) {
        op.lap(i) = 0;
    }
    for (int i : domain.boundary()) {
        op.value(i) = 1;
    }
    /// [Implicit syntax]

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
    solver.compute(M);
    Eigen::VectorXd sol = solver.solve(rhs);
    for (int i = 0; i < domain.size(); ++i) {
        EXPECT_NEAR(sol[i], 1.0, 1e-14);
    }
}

TEST(Operators, Implicit2dGradOfScalar) {
    /// [Implicit Grad example]
    // Solve problem for S unknown scalar field:
    // grad S . v + laplace S = phi
    // grad S . (y, x) + alfa laplace S = 4 (alfa + xy) + e^x(2 x cos(2y) - (3 alfa - y) sin(2y))
    // on [0, 1] x [0, 1] with BC:
    // S(x, 0) = x^2
    // S(x, 1) = 1 + x^2 + sin(2) e^x
    // S(0, y) = y^2 + sin(2y)
    // S(1, y) = 1 + y^2 + e * sin(2y)
    // with solution:
    // S(x, y) = e^x sin(2y) + x^2 + y^2
    std::function<Vec2d(Vec2d)> analytical = [] (const Vec2d& p) {
        return std::exp(p[0]) * std::sin(2*p[1]) + p[0]*p[0] + p[1]*p[1];
    };
    // Prepare domain
    BoxShape<Vec2d> box(0.0, 1.0);
    DomainDiscretization<Vec2d> domain = box.discretizeWithStep(0.05);
    int N = domain.size();
    domain.findSupport(FindClosest(9));
    // Prepare operators and matrix
    WLS<Monomials<Vec2d>, GaussianWeight<Vec2d>> mls(2, 15.0 / N);
    auto shapes = domain.computeShapes<sh::lap|sh::grad>(mls, domain.interior());
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    M.reserve(shapes.supportSizes());
    Eigen::VectorXd rhs(N); rhs.setZero();

    auto op = shapes.implicitOperators(M, rhs);

    double alpha = 0.25;
    for (int i : domain.interior()) {
        double x = domain.pos(i, 0), y = domain.pos(i, 1);
        op.grad(i, {y, x}) + alpha * op.lap(i) =
                4*(alpha + x*y) + std::exp(x) * (2*x*std::cos(2*y) - (3*alpha - y)*std::sin(2*y));
    }
    // Set boundary conditions
    for (int i : domain.types() == -1) {
        double y = domain.pos(i, 1);
        op.value(i) = y*y + std::sin(2*y);
    }
    for (int i : domain.types() == -2) {
        double y = domain.pos(i, 1);
        op.value(i) = 1 + y*y + std::exp(1) * std::sin(2*y);
    }
    for (int i : domain.types() == -3) {
        double x = domain.pos(i, 0);
        op.value(i) = x*x;
    }
    for (int i : domain.types() == -4) {
        double x = domain.pos(i, 0);
        op.value(i) = 1 + x*x + std::sin(2) * std::exp(x);
    }

    // Sparse solve
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
    solver.compute(M);
    Eigen::VectorXd sol = solver.solve(rhs);
    for (int i = 0; i < N; ++i) {
        Vec2d correct = analytical(domain.pos(i));
        ASSERT_NEAR(sol[i], correct[0], 5e-4);
    }
    /// [Implicit Grad example]
}

TEST(Operators, ImplicitDer1) {
    /// [Implicit derivative example]
    // Solve problem for F unknown vector field:
    // f''(x) = 6x, f'(1) = 1, f(0) = 2
    // with solution:
    // f(x) = 2 - 2x + x^3
    std::function<Vec1d(Vec1d)> analytical = [] (const Vec1d& p) {
        return Vec1d(2 - 2*p[0] + p[0]*p[0]*p[0]);
    };
    // Prepare domain
    BoxShape<Vec1d> box(0.0, 1.0);
    double dx = 0.02;
    DomainDiscretization<Vec1d> domain = box.discretizeWithStep(dx);
    int N = domain.size();
    int support_size = 3;
    domain.findSupport(FindClosest(support_size));
    // Prepare operators and matrix
    WLS<Monomials<Vec1d>, GaussianWeight<Vec1d>> wls(2, 15.0 / N);

    auto storage = domain.computeShapes<sh::d1|sh::lap>(wls);

    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(N, N);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(N);

    ImplicitOperators<decltype(storage), decltype(M), decltype(rhs)> op(storage, M, rhs);

    // Set equation on interior
    for (int i : domain.interior()) {
        op.lap(i) = 6*domain.pos(i, 0);
    }
    // Set boundary conditions
    for (int i : domain.types() == -1) {
        op.value(i) = 2;
    }
    for (int i : domain.types() == -2) {
        op.der1(i, 0) = 1;
    }

    Eigen::VectorXd sol = M.lu().solve(rhs);
    for (int i = 0; i < N; ++i) {
        Vec1d correct = analytical(domain.pos(i));
        ASSERT_NEAR(sol[i], correct[0], 1e-3);
    }
    /// [Implicit derivative example]
}

TEST(Operators, DISABLED_ImplicitOperatorsUsageExample) {
    // Prepare domain
    BoxShape<Vec2d> box(0.0, 1.0);
    DomainDiscretization<Vec2d> domain = box.discretizeWithStep(0.05);
    int N = domain.size();
    domain.findSupport(FindClosest(9));
    // Prepare operators and matrix
    WLS<Monomials<Vec2d>, GaussianWeight<Vec2d>> mls(2, 15.0 / N);
    auto shapes = domain.computeShapes(mls, domain.interior());

    /// [Implicit operators usage example]
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    M.reserve(shapes.supportSizes());
    Eigen::VectorXd rhs(N); rhs.setZero();

    auto op = shapes.implicitOperators(M, rhs);
    std::cout << op << std::endl;

    double alpha = 0.25;
    for (int i : domain.interior()) {
        double x = domain.pos(i, 0), y = domain.pos(i, 1);
        op.grad(i, {y, x}) + alpha * op.lap(i) =
                4*(alpha + x*y) + std::exp(x) * (2*x*std::cos(2*y) - (3*alpha - y)*std::sin(2*y));
    }
    // Set boundary conditions
    for (int i : domain.types() == -1) {
        double y = domain.pos(i, 1);
        op.value(i) = y*y + std::sin(2*y);
    }
    for (int i : domain.types() == -2) {
        double y = domain.pos(i, 1);
        op.value(i) = 1 + y*y + std::exp(1) * std::sin(2*y);
    }
    for (int i : domain.types() == -3) {
        double x = domain.pos(i, 0);
        op.value(i) = x*x;
    }
    for (int i : domain.types() == -4) {
        double x = domain.pos(i, 0);
        op.value(i) = 1 + x*x + std::sin(2) * std::exp(x);
    }
    /// [Implicit operators usage example]
}

}  // namespace mm
