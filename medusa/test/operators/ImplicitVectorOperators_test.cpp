#include <medusa/bits/approximations/Gaussian_fwd.hpp>
#include <medusa/bits/approximations/JacobiSVDWrapper_fwd.hpp>
#include <medusa/bits/approximations/Monomials.hpp>
#include <medusa/bits/approximations/WeightFunction_fwd.hpp>
#include <medusa/bits/approximations/WLS_fwd.hpp>
#include <medusa/bits/domains/BoxShape_fwd.hpp>
#include <medusa/bits/domains/DomainDiscretization_fwd.hpp>
#include <medusa/bits/domains/FindClosest.hpp>
#include <medusa/bits/operators/computeShapes.hpp>
#include <medusa/bits/operators/ImplicitVectorOperators.hpp>
#include <medusa/bits/operators/RaggedShapeStorage.hpp>
#include <medusa/bits/operators/ShapeStorage.hpp>
#include <medusa/bits/operators/UniformShapeStorage.hpp>
#include <medusa/bits/types/Vec.hpp>
#include <cmath>

#include "gtest/gtest.h"
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

namespace mm {

class ImplicitVecOp : public ::testing::Test {
  public:
    UniformShapeStorage<Vec2d> shapes;
    Eigen::Matrix<double, 6, 6> M;
    Eigen::Matrix<double, 6, 1> rhs;
    ImplicitVectorOperators<decltype(shapes), decltype(M), decltype(rhs)> op;

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

        op = shapes.implicitVectorOperators(M, rhs);
    }
};

TEST_F(ImplicitVecOp, AssignTwice) {
    auto expectedM = M;
    auto expectedRhs = rhs;

    expectedM(0, 0) = 1;
    expectedM(3, 3) = 1;
    expectedRhs[0] = 1;
    expectedRhs[3] = 1;
    op.value(0) = 1;
    EXPECT_EQ(0, (expectedM - M).lpNorm<1>());
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());

    expectedM(0, 0) = 2;
    expectedM(3, 3) = 2;
    expectedRhs[0] = 2;
    expectedRhs[3] = 2;
    op.value(0) = 1;
    EXPECT_EQ(0, (expectedM - M).lpNorm<1>());
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());
}

TEST_F(ImplicitVecOp, AssignCombined) {
    double v = 4.2;
    auto expectedM = M;
    auto expectedRhs = rhs;
    op.value(1) + v * op.lap(1) + op.value(1) = {v, 2*v};
    expectedM.row(1).head(3) += v * shapes.laplace(1);
    expectedM.row(4).tail(3) += v * shapes.laplace(1);
    expectedM(1, 1) += 2;
    expectedM(4, 4) += 2;
    expectedRhs[1] = v;
    expectedRhs[4] = 2*v;
    EXPECT_EQ(0, (expectedM - M).lpNorm<1>());
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());
}

TEST_F(ImplicitVecOp, AssignValue) {
    double v = 4.2, f = -3.4;
    auto expectedM = M;
    auto expectedRhs = rhs;
    expectedM(1, 1) = f;
    expectedM(4, 4) = f;
    expectedRhs[1] = v;
    expectedRhs[4] = 2*v;
    f*op.value(1) = {v, 2*v};
    EXPECT_EQ(0, (expectedM - M).lpNorm<1>());
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());
}

TEST_F(ImplicitVecOp, AssignLap) {
    double v = 4.2;
    auto expectedM = M;
    auto expectedRhs = rhs;
    expectedM.row(0).head(3) = shapes.laplace(0);
    expectedM.row(3).tail(3) = shapes.laplace(0);
    expectedRhs[0] = v;
    expectedRhs[3] = 2*v;
    op.lap(0) = {v, 2*v};
    EXPECT_EQ(0, (expectedM - M).lpNorm<1>());
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());
}

TEST_F(ImplicitVecOp, AssignGrad) {
    Vec2d direction(3.4, 2.4);
    double v = 4.2, v_row = 5.3;
    auto expectedM = M;
    auto expectedRhs = rhs;

    expectedM.row(0).head(3) = direction[0]*shapes.d1(0, 0) + direction[1]*shapes.d1(1, 0);
    expectedM.row(1).head(3) = direction[0]*shapes.d1(0, 0) + direction[1]*shapes.d1(1, 0);
    expectedM.row(3).tail(3) = direction[0]*shapes.d1(0, 0) + direction[1]*shapes.d1(1, 0);
    expectedM.row(4).tail(3) = direction[0]*shapes.d1(0, 0) + direction[1]*shapes.d1(1, 0);
    expectedRhs[0] = v;
    expectedRhs[1] = v_row;
    expectedRhs[3] = 2*v;
    expectedRhs[4] = 2*v_row;
    op.grad(0, direction) = {v, 2*v};
    op.grad(0, direction, 1) = {v_row, 2*v_row};
    EXPECT_EQ(0, (expectedM - M).lpNorm<1>());
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());
}

TEST_F(ImplicitVecOp, AssignNeumann) {
    Vec2d direction(0.6, 0.8);
    double v = 4.2, v_row = 5.3;
    auto expectedM = M;
    auto expectedRhs = rhs;

    expectedM.row(0).head(3) = direction[0]*shapes.d1(0, 0) + direction[1]*shapes.d1(1, 0);
    expectedM.row(1).head(3) = direction[0]*shapes.d1(0, 0) + direction[1]*shapes.d1(1, 0);
    expectedM.row(3).tail(3) = direction[0]*shapes.d1(0, 0) + direction[1]*shapes.d1(1, 0);
    expectedM.row(4).tail(3) = direction[0]*shapes.d1(0, 0) + direction[1]*shapes.d1(1, 0);
    expectedRhs[0] = v;
    expectedRhs[1] = v_row;
    expectedRhs[3] = 2*v;
    expectedRhs[4] = 2*v_row;
    op.neumann(0, direction) = {v, 2*v};
    op.neumann(0, direction, 1) = {v_row, 2*v_row};
    EXPECT_EQ(0, (expectedM - M).lpNorm<1>());
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());
}

TEST_F(ImplicitVecOp, AssignGradDiv) {
    double v = 4.2;
    auto expectedM = M;
    auto expectedRhs = rhs;

    expectedM.row(0).head(3) = shapes.d2(0, 0, 0);  // dxdx in node 0
    expectedM.row(0).tail(3) = shapes.d2(0, 1, 0);  // dxdx in node 0
    expectedM.row(3).head(3) = shapes.d2(0, 1, 0);  // dxdy in node 1
    expectedM.row(3).tail(3) = shapes.d2(1, 1, 0);  // dydy in node 2
    expectedRhs[0] = v;
    expectedRhs[3] = v;
    op.graddiv(0) = v;

    EXPECT_EQ(0, (expectedM - M).lpNorm<1>());
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());
}

TEST_F(ImplicitVecOp, AssignDeath) {
    EXPECT_DEATH(op.lap(0) + op.lap(1) = 5, "Cannot add together terms for different matrix rows");
}

TEST_F(ImplicitVecOp, Syntax) {
    0.0*op.lap(0, 0) + op.lap(0, 0) = 1;
    op.lap(0, 1) + 0.0*op.lap(0, 1) = 1;
    0.25*op.lap(0, 2) + 0.75*op.lap(0, 2) = 1;

    for (int i = 1; i < 3; ++i) {
        EXPECT_EQ(M.row(i), M.row(0));
        EXPECT_EQ(rhs.row(i), rhs.row(0));
    }
}


TEST(Operators, ImplicitVectorComplexNumbers) {
    UniformShapeStorage<Vec2d> shapes;
    shapes.resize({2, 2});

    shapes.setSupport(0, {0, 1});
    shapes.setSupport(1, {0, 1});

    Eigen::Vector2d s0(-0.12, 2.42), s1(0.34, 23.3);
    std::vector<Eigen::Vector2d> s = {s0, s1};
    for (int i = 0; i < 2; ++i) {
        shapes.setLaplace(i, s[i]);
        for (int d1 = 0; d1 < 2; ++d1) {
            shapes.setD1(d1, i, s[i]);
            for (int d2 = d1; d2 < 2; ++d2) {
                shapes.setD2(d1, d2, i, s[i]);
            }
        }
    }

    Eigen::Matrix4cd M; M.setZero();
    Eigen::Vector4cd rhs; rhs.setZero();

    auto op = shapes.implicitVectorOperators(M, rhs);
    std::complex<double> c1(1.2, -3.4), c2(0.7, 8.9), c3(-5.6, 2.3);
    double x = 1.6;
    Vec2d dir(3.4, 2.3);
    c1*op.lap(0) + c2*op.value(0) = {c3, c1};
    op.grad(1, dir) + c1*op.graddiv(1) + x*op.lap(1) = {c3, c1};

    // eq 1
    Eigen::Matrix4cd expectedM = M;
    Eigen::Vector4cd expectedRhs = rhs;
    expectedM.row(0).head(2) = c1*s0;
    expectedM(0, 0) += c2;
    expectedRhs[0] = c3;
    expectedM.row(2).tail(2) = c1*s0;
    expectedM(2, 2) += c2;
    expectedRhs[2] = c1;

    EXPECT_LT((expectedM - M).lpNorm<1>(), 3e-16);
    EXPECT_EQ(0, (expectedRhs - rhs).lpNorm<1>());
}

TEST(Operators, implicit2dLaplaceOfvector) {
    /// [Implicit vector laplace 2d example]
    // Solve problem for F unknown vector field:
    // (Laplace F)(x, y) = [4 (1 - 2 pi^2 y^2) sin(2 pi x); 4 (1 - x - x y) e^(2 y)]
    // on [0, 1] x [0, 1] with BC:
    // F(x, 0) = [0, 1]
    // F(x, 1) = [2 sin(2 pi x), (1 - x) e^2]
    // F(0, y) = [0, e^(2 y)]
    // F(1, y) = [0, (1 - y) e^(2 y)]
    // with solution:
    // F(x, y) = [2 sin(2 pi x) y^2,  (1 - x y) e^(2 y)]
    std::function<Vec2d(Vec2d)> analytical = [] (const Vec2d& p) {
        return Vec2d({2*std::sin(2*PI*p[0])*p[1]*p[1],  (1 - p[0]*p[1]) * std::exp(2*p[1])});
    };
    // Prepare domain
    BoxShape<Vec2d> box(0., 1.);
    DomainDiscretization<Vec2d> domain = box.discretizeWithStep(0.025);
    int N = domain.size();
    int support_size = 9;
    domain.findSupport(FindClosest(support_size));
    // Prepare operators and matrix
    WLS<Monomials<Vec2d>, GaussianWeight<Vec2d>> wls(2, std::sqrt(15.0 / N));
    auto storage = domain.computeShapes(wls);
    Eigen::SparseMatrix<double> M(2*N, 2*N);
    M.reserve(storage.supportSizes()+storage.supportSizes());
    Eigen::VectorXd rhs(2*N); rhs.setZero();
    auto op = storage.implicitVectorOperators(M, rhs);
    // Set equation on interior
    for (int i : domain.interior()) {
        double x = domain.pos(i, 0), y = domain.pos(i, 1);
        op.lap(i) = {4 * (1 - 2 * PI*PI * y*y) * std::sin(2 * PI * x),
                     4 * (1 - x - x * y) * std::exp(2 * y)};
    }
    // Set boundary conditions
    for (int i : domain.types() == -1) {
        double y = domain.pos(i, 1);
        op.value(i) = {0, std::exp(2*y)};
    }
    for (int i : domain.types() == -2) {
        double y = domain.pos(i, 1);
        op.value(i) = {0, (1-y) * std::exp(2*y)};
    }
    for (int i : domain.types() == -3) {
        op.value(i) = {0, 1};
    }
    for (int i : domain.types() == -4) {
        double x = domain.pos(i, 0);
        op.value(i) = {2 * std::sin(2 * PI * x), (1 - x) * std::exp(2)};
    }
    // Sparse solve
    M.makeCompressed();

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(M);
    Eigen::VectorXd sol = solver.solve(rhs);
    for (int i = 0; i < N; ++i) {
        Vec2d correct = analytical(domain.pos(i));
        EXPECT_NEAR(sol[i], correct[0], 1e-3);
        EXPECT_NEAR(sol[N+i], correct[1], 1e-3);
    }
    /// [Implicit vector laplace 2d example]
}

TEST(Operators, Implicit2dGradOfVector) {
    /// [Implicit Gradvec example]
    // Solve problem for u unknown vector field:
    // grad u . v + 1/2 laplace u = f
    // grad u . (y, x) + 1/2 laplace u = (x^3+y+2 x y^2, x^2+(y-1) y)
    // on [0, 1] x [0, 1] with BC:
    // u(x, 0) = (0, -x)
    // u(x, 1) = (x*x, 0)
    // u(0, y) = (0, 0)
    // u(1, y) = (y, y-1)
    // with solution:
    // u(x, y) = (x^2 y, yx - x)
    std::function<Vec2d(Vec2d)> analytical = [] (const Vec2d& p) {
        return Vec2d(p[0]*p[0]*p[1], (p[1]-1)*p[0]);
    };
    // Prepare domain
    BoxShape<Vec2d> box(0., 1.);
    DomainDiscretization<Vec2d> domain = box.discretizeWithStep(0.05);
    int N = domain.size();
    int support_size = 9;
    domain.findSupport(FindClosest(support_size));
    // Prepare operators and matrix
    WLS<Monomials<Vec2d>, GaussianWeight<Vec2d>> wls(2, 15.0 / N);
    auto storage = domain.computeShapes(wls);
    Eigen::SparseMatrix<double> M(2*N, 2*N);
    M.reserve(storage.supportSizes()+storage.supportSizes());
    Eigen::VectorXd rhs(2*N); rhs.setZero();
    auto op = storage.implicitVectorOperators(M, rhs);

    Range<Vec2d> v(N);
    for (int i = 0; i < N; ++i) {
        v[i] = Vec2d({domain.pos(i, 1), domain.pos(i, 0)});
    }
    // Set equation on interior
    for (int i : domain.interior()) {
        double x = domain.pos(i, 0), y = domain.pos(i, 1);
        op.grad(i, v[i]) + 0.5*op.lap(i) = {x*x*x + y + 2*x*y*y, x*x + (y-1)*y};
    }
    // Set boundary conditions
    for (int i : domain.types() == -1) {
        op.value(i) = 0;
    }
    for (int i : domain.types() == -2) {
        double y = domain.pos(i, 1);
        op.value(i) = {y, y-1};
    }
    for (int i : domain.types() == -3) {
        double x = domain.pos(i, 0);
        op.value(i) = {0, -x};
    }
    for (int i : domain.types() == -4) {
        double x = domain.pos(i, 0);
        op.value(i) = {x*x, 0};
    }

    // Sparse solve
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
    solver.compute(M);
    Eigen::VectorXd sol = solver.solve(rhs);
    for (int i = 0; i < N; ++i) {
        Vec2d correct = analytical(domain.pos(i));
        ASSERT_NEAR(sol[i], correct[0], 1e-3);
        ASSERT_NEAR(sol[i+N], correct[1], 1e-3);
    }
    /// [Implicit Gradvec example]
}

TEST(Operators, Implicit2dGradDiv) {
    /// [Implicit graddiv example]
    // Solve problem for u unknown vector field:
    // grad div u + laplace u = (2 + 4xy - 8pi^2sin(2*pi*x), 2(2x^2+y^2))
    // on [0, 1] x [0, 1] with BC:
    // u(0, y) = (y^2, 0)
    // u(1, y) = (y^2, y^2)
    // u(x, 0) = (sin(2*pi*x), 0)
    // u(x, 1) = (1+sin(2*pi*x, x^2)
    // with solution:
    // u(x, y) = (y^2 + sin(2*pi*x), x^2y^2)
    std::function<Vec2d(Vec2d)> analytical = [] (const Vec2d& p) {
        return Vec2d({p[1]*p[1]+std::sin(2*PI*p[0]), p[0]*p[0]*p[1]*p[1]});
    };
    // Prepare domain
    BoxShape<Vec2d> box(0., 1.);
    DomainDiscretization<Vec2d> domain = box.discretizeWithStep(0.025);
    int N = domain.size();
    int support_size = 6;
    domain.findSupport(FindClosest(support_size));
    // Prepare operators and matrix
    WLS<Monomials<Vec2d>, GaussianWeight<Vec2d>> wls(2, std::sqrt(15.0 / N));
    auto storage = domain.computeShapes<sh::graddiv|sh::lap>(wls, domain.interior());

    Eigen::SparseMatrix<double, Eigen::RowMajor> M(2*N, 2*N);
    M.reserve(Range<int>(2*N, support_size));
    Eigen::VectorXd rhs(2*N); rhs.setZero();

    auto op = storage.implicitVectorOperators(M, rhs);
    // Set equation on interior
    for (int i : domain.interior()) {
        double x = domain.pos(i, 0), y = domain.pos(i, 1);
        op.lap(i) + op.graddiv(i) = {2+4*x*y-8*PI*PI*std::sin(2*PI*x), 2*(2*x*x + y*y)};
    }
    // Set boundary conditions
    for (int i : domain.types() == -1) {
        double y = domain.pos(i, 1);
        op.value(i) = {y*y, 0};
    }
    for (int i : domain.types() == -2) {
        double y = domain.pos(i, 1);
        op.value(i) = {y*y, y*y};
    }
    for (int i : domain.types() == -3) {
        double x = domain.pos(i, 0);
        op.value(i) = {std::sin(2*PI*x), 0};
    }
    for (int i : domain.types() == -4) {
        double x = domain.pos(i, 0);
        op.value(i) = {1+std::sin(2*PI*x), x*x};
    }

    // Sparse solve
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
    solver.compute(M);
    Eigen::VectorXd sol = solver.solve(rhs);
    for (int i = 0; i < N; ++i) {
        Vec2d correct = analytical(domain.pos(i));
        EXPECT_NEAR(sol[i], correct[0], 0.8e-2);
        EXPECT_NEAR(sol[i+N], correct[1], 0.8e-2);
    }
    /// [Implicit graddiv example]
}

TEST(Operators, ImplicitDer1Sum) {
    /// [Sum of derivatives example]
    // The vector field is given by:
    // u = u_0*(sin(pi*(x^2+y^2))+1)
    // v = v_0*(cos(pi*(x^2+y^2))+1)
    // The derivatives are given by:
    // u_x = 2*pi*u_0*x*cos(pi*(x^2+y^2))
    // u_y = 2*pi*u_0*y*cos(pi*(x^2+y^2))
    // v_x = -2*pi*v_0*x*sin(pi*(x^2+y^2))
    // v_y = -2*pi*v_0*y*sin(pi*(x^2+y^2))

    double u0 = 1, v0 = 1;
    double alpha = 1, beta = 1, gamma = 1, delta = 1;

    // Vector field function that returns (u,v)
    std::function<Vec2d(Vec2d)> analytical = [&u0, &v0] (const Vec2d& p) {
        return Vec2d({u0*(std::sin(p.squaredNorm())+1),
                      v0*(std::cos(p.squaredNorm())+1)});
    };
    // Analytical function for \alpha*u_x + \beta*v_y
    std::function<double(Vec2d)> tx = [&alpha, &beta, &u0, &v0] (const Vec2d& p) {
        return alpha*(2*u0*p[0]*std::cos(p.squaredNorm())) +
               beta*(-2*v0*p[1]*std::sin(p.squaredNorm()));
    };
    // Analytical function for \gamma*u_y + \delta*v_x
    std::function<double(Vec2d)> ty = [&gamma, &delta, &u0, &v0] (const Vec2d& p) {
        return gamma*(2*u0*p[1]*std::cos(p.squaredNorm())) +
               delta*(-2*v0*p[0]*std::sin(p.squaredNorm()));
    };
    // Create and fill domain
    BoxShape<Vec2d> box({-0.2, -0.5}, {0.8, 0.5});
    double step = 0.01;
    DomainDiscretization<Vec2d> domain = box.discretizeWithStep(step);
    int support_size = 12;
    domain.findSupport(FindClosest(support_size));
    int N = domain.size();
    // Prepare MLS engine
    WLS<Monomials<Vec2d>, GaussianWeight<Vec2d>> wls(2, step);
    auto storage = domain.computeShapes<sh::d1>(wls);
    // Fill Matrix with operators for
    // \alpha*u_x + \beta*v_y in x-direction and
    // \gamma*u_y + \delta*v_x in y-direction
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(2 * N, 2 * N);
    M.reserve(Range<int>(2 * N, 2 * support_size));
    Eigen::VectorXd field(2 * N), rhs(2 * N);
    rhs.setZero();
    auto op = storage.implicitVectorOperators(M, rhs);

    for (int i = 0; i < N; ++i) {
        // apply in x-direction
        /// [Eq usage example]
        alpha*op.eq(0).c(0).der1(i, 0) + beta*op.eq(0).c(1).der1(i, 1) = tx(domain.pos(i));
        /// [Eq usage example]
        // apply in y-direction
        gamma*op.eq(1).c(0).der1(i, 1) + delta*op.eq(1).c(1).der1(i, 0) = ty(domain.pos(i));

        auto f = analytical(domain.pos(i));
        field(i) = f[0];
        field(i+N) = f[1];
    }
    // perform multiplication
    Eigen::VectorXd prod = M*field;

    // check both components are close enough
    for (int i = 0; i < N; ++i) {
        ASSERT_NEAR(prod[i], rhs[i], 1e-3);
        ASSERT_NEAR(prod[i+N], rhs[i+N], 1e-3);
    }
    /// [Sum of derivatives example]
}

TEST(Operators, DISABLED_ImplicitVectorOperatorsUsageExample) {
    BoxShape<Vec2d> box({-0.2, -0.5}, {0.8, 0.5});
    double step = 0.01;
    DomainDiscretization<Vec2d> domain = box.discretizeWithStep(step);
    int support_size = 12;
    domain.findSupport(FindClosest(support_size));
    int N = domain.size();
    // Prepare MLS engine
    WLS<Monomials<Vec2d>, GaussianWeight<Vec2d>> wls(2, step);
    auto storage = domain.computeShapes(wls);

    /// [Implicit vector operators usage example]
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(2 * N, 2 * N);
    M.reserve(Range<int>(2 * N, 2 * support_size));
    Eigen::VectorXd field(2 * N), rhs(2 * N);

    auto op = storage.implicitVectorOperators(M, rhs);
    std::cout << op << std::endl;

    // Set equation on interior
    for (int i : domain.interior()) {
        double x = domain.pos(i, 0), y = domain.pos(i, 1);
        op.lap(i) + op.graddiv(i) = {2+4*x*y-8*PI*PI*std::sin(2*PI*x), 2*(2*x*x + y*y)};
    }
    // Set boundary conditions
    for (int i : domain.types() == -1) {
        double y = domain.pos(i, 1);
        op.value(i) = {y*y, 0};
    }
    for (int i : domain.types() == -2) {
        double y = domain.pos(i, 1);
        op.value(i) = {y*y, y*y};
    }
    for (int i : domain.types() == -3) {
        double x = domain.pos(i, 0);
        op.value(i) = {std::sin(2*PI*x), 0};
    }
    for (int i : domain.types() == -4) {
        double x = domain.pos(i, 0);
        op.value(i) = {1+std::sin(2*PI*x), x*x};
    }
    /// [Implicit vector operators usage example]
}

}  // namespace mm
