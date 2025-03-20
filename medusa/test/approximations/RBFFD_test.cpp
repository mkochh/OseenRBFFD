#include <medusa/bits/approximations/RBFFD.hpp>
#include <medusa/bits/approximations/Gaussian.hpp>
#include <medusa/bits/approximations/Operators.hpp>
#include <medusa/bits/approximations/Monomials.hpp>
#include <medusa/bits/approximations/Polyharmonic.hpp>
#include <medusa/bits/approximations/RBFBasis.hpp>
#include <medusa/bits/types/Vec.hpp>
#include <Eigen/Cholesky>
#include <Eigen/LU>

#include "gtest/gtest.h"

namespace mm {

TEST(Approximations, RBFFDGauss2D) {
    double h = 0.1123;
    double s = 1.5;

    Gaussian<double> g(s);
    RBFFD<Gaussian<double>, Vec2d> appr(g);

    Range<Vec2d> support = {{0, 0}, {0, h}, {h, 0}, {0, -h}, {-h, 0},
                            {-h, h}, {-h, -h}, {h, -h}, {h, h}};
    appr.compute({0.0, 0.0}, support);

    double a = -4.*(s*s + h*h/std::pow(std::sinh((h/s)*(h/s)), 2)) / std::pow(s, 4);
    double b = 4.*std::exp(3.*(h/s)*(h/s))*h*h / std::pow(-1+std::exp(2*(h/s)*(h/s)), 2)
               / std::pow(s, 4);

    Eigen::VectorXd shape = appr.getShape(Lap<2>());
    Eigen::VectorXd expected(9); expected << a, b, b, b, b, 0, 0, 0, 0;
    ASSERT_EQ(expected.size(), shape.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_NEAR(expected[i], shape[i], 5e-6);
    }

    double c = 2.*std::exp(3.*(h/s)*(h/s))*h / (-1+std::exp(4*(h/s)*(h/s))) / s / s;
    shape = appr.getShape(Der1<2>(1));
    expected << 0, c, 0, -c, 0, 0, 0, 0, 0;
    ASSERT_EQ(expected.size(), shape.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_NEAR(expected[i], shape[i], 5e-8);
    }
}

TEST(Approximations, RBFFDGauss2DAugConst) {
    double h = 0.1123;
    double s = 1.5;

    Gaussian<double> g(s);
    RBFFD<Gaussian<double>, Vec2d> appr(g, 0);

    Range<Vec2d> support = {{0, 0}, {0, h}, {h, 0}, {0, -h}, {-h, 0},
                            {-h, h}, {-h, -h}, {h, -h}, {h, h}};
    appr.compute({0.0, 0.0}, support);

    Eigen::VectorXd shape = appr.getShape(Lap<2>());
    // computed with Mathematica from (quite long) analytical expression
    Eigen::VectorXd expected(9);
    expected << -336.1792004798, 88.49955236455, 88.49955236455, 88.49955236455, 88.49955236455,
            -4.4547522446, -4.4547522446, -4.4547522446, -4.4547522446;
    ASSERT_EQ(expected.size(), shape.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_NEAR(expected[i], shape[i], 5e-6);
    }

    shape = appr.getShape(Der1<2>(1));
    // same as without a constant
    double c = 2.*std::exp(3.*(h/s)*(h/s))*h / (-1+std::exp(4*(h/s)*(h/s))) / s / s;
    expected << 0, c, 0, -c, 0, 0, 0, 0, 0;
    ASSERT_EQ(expected.size(), shape.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_NEAR(expected[i], shape[i], 5e-8);
    }
}

TEST(Approximations, RBFFDPhs) {
    Polyharmonic<double, 3> phs;
    Monomials<Vec2d> mon(2);
    RBFFD<decltype(phs), Vec2d, NoScale, Eigen::PartialPivLU<Eigen::MatrixXd>> appr(phs, mon);
    double h = 0.1234;
    Range<Vec2d> supp = {{0, 0}, {0, h}, {h, 0}, {0, -h}, {-h, 0},
                         {-h, h}, {-h, -h}, {h, -h}, {h, h}};

    double a = -((356 + 213*std::sqrt(2) + 40*std::sqrt(5) + 275*std::sqrt(10))/1312.0);
    Eigen::VectorXd expected(supp.size());
    expected << 4 * a - 4, -2 * a + 1, -2 * a + 1, -2 * a + 1, -2 * a + 1, a, a, a, a;
    expected /= h * h;

    appr.compute({0.0, 0.0}, supp);
    auto sh = appr.getShape(Lap<2>());
    for (int i = 0; i < supp.size(); ++i) {
        EXPECT_NEAR(expected[i], sh[i], 1e-11);
    }

    /**
     * Test computed with Mathematica:
     *
     * phig = (Exp[-#/sigma^2] &);
     * phip = (#^(3/2) &);
     *
     * $Assumptions = {h > 0, sigma > 0, c > 0, x0 < x, phi[-h] == phi[h],
     *    phi'[-h] == phi'[h], phi''[-h] == phi''[h]};
     *
     * pts = {{0, 0}, {0, h}, {h, 0}, {0, -h}, {-h, 0}, {-h,
     *     h}, {-h, -h}, {h, -h}, {h, h}};
     * c = {0, 0};
     * phir[x_] := phi[(x - c).(x - c)/h^2];
     * n = Length[pts];
     *
     * A = Table[phir[pts[[i]] - pts[[j]]], {i, 1, n}, {j, 1, n}]  //
     *    Simplify;
     * p = {(1 &), (#[[2]] &), (#[[2]]^2 &),  (#[[1]] &), (#[[1]] #[[2]] &),
     *      (#[[1]]^2 &)}; (* , (#^3&), (#^4&)*)
     * P = Table[f[(x - c)/h], {x, pts}, {f, p}] ;
     * s = Length[p];
     * M = Table[0, {i, 1, n + s}, {j, 1, n + s}];
     * M[[1 ;; n, 1 ;; n]] = A;
     * M[[n + 1 ;; n + s, 1 ;; n]] = Transpose[P];
     * M[[1 ;; n, n + 1 ;; n + s]] = P;
     * M // MatrixForm
     *
     * coef = Table[al[i], {i, n}]~Join~Table[be[i], {i, s}];
     * rhsphi = Table[
     *      Laplacian[phir[{x, y} - pts[[i]]], {x, y}], {i, 1, n}] /. {x -> c[[1]], y -> c[[2]]};
     * rhsmon = Table[
     *      Laplacian[f[({x, y} - c)/h], {x, y}], {f, p}] /. {x -> c[[1]], y -> c[[2]]};
     * rhs = Join[rhsphi, rhsmon] // Simplify;
     *
     * Mp = M /. phi -> phip // Simplify;
     * rhsp = rhs /. phi -> phip  // Simplify;
     * MatrixForm /@ {Mp, rhsp}
     *
     * sol = LinearSolve[Mp, rhsp] // FullSimplify;
     * cf = sol[[1 ;; n]];
     * uu = Table[u[p], {p, pts}];
     * a = -((356 + 213 Sqrt[2] + 40 Sqrt[5] + 275 Sqrt[10])/1312);
     * cf2 = 1/h^2 {4 a - 4, -2 a + 1, -2 a + 1, -2 a + 1, -2 a + 1, a, a, a, a};
     * cf2 - cf // Simplify
     */
}

TEST(Approximations, RBFFDGauss2DAugOrd1) {
    double h = 0.1123;
    double s = 1.5;

    Gaussian<double> g(s);
    RBFFD<Gaussian<double>, Vec2d> appr(g, 1);

    Range<Vec2d> support = {{0, 0}, {0, h}, {h, 0}, {0, -h}, {-h, 0},
                            {-h, h}, {-h, -h}, {h, -h}, {h, h}};
    appr.compute({0.0, 0.0}, support);

    Eigen::VectorXd shape = appr.getShape(Lap<2>());
    // computed with Mathematica from (quite long) analytical expression (same as before)
    Eigen::VectorXd expected(9);
    expected << -336.1792004798, 88.49955236455, 88.49955236455, 88.49955236455, 88.49955236455,
            -4.4547522446, -4.4547522446, -4.4547522446, -4.4547522446;
    ASSERT_EQ(expected.size(), shape.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_NEAR(expected[i], shape[i], 5e-6);
    }

    shape = appr.getShape(Der1<2>(1));
    // computed with Mathematica from (quite long) analytical expression
    expected << 0, 5.94478246428, 0, -5.94478246428, 0, -0.74621135681,
            0.74621135681, 0.74621135681, -0.74621135681;
    ASSERT_EQ(expected.size(), shape.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_NEAR(expected[i], shape[i], 5e-8);
    }
}

TEST(Approximations, DISABLED_RBFFDUsageExmaple) {
    /// [RBFFD usage example]
    double h = 0.1123;
    double s = 1.5;

    Gaussian<double> g(s);
    // Gaussian RBF's augmented with a constant
    RBFFD<Gaussian<double>, Vec2d> appr(g, Monomials<Vec2d>(0));
    std::cout << appr << std::endl;

    // Local neighbourhood
    Range<Vec2d> support = {{0, 0}, {0, h}, {h, 0}, {0, -h}, {-h, 0},
                            {-h, h}, {-h, -h}, {h, -h}, {h, h}};

    // compute approximations at point `{0.0, 0.0}`.
    appr.compute({0.0, 0.0}, support);

    // Get info about the computation.
    std::cout << appr.rbf() << std::endl;
    std::cout << appr.monomials() << std::endl;
    std::cout << appr.center() << std::endl;
    std::cout << appr.scale() << std::endl;
    std::cout << appr.localCoordinates() << std::endl;
    std::cout << appr.getMatrix() << std::endl;
    Eigen::PartialPivLU<Eigen::MatrixXd> solver = appr.solver();


    // Get shape (stencil weights) for approximation of Laplacian.
    Eigen::VectorXd shape = appr.getShape(Lap<2>());

    // Get shape (stencil weights) for approximation of value
    shape = appr.getShape();
    /// [RBFFD usage example]
    (void) solver;
}

}  // namespace mm
