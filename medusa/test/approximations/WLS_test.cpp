#include <medusa/bits/approximations/WLS.hpp>
#include <medusa/bits/approximations/JacobiSVDWrapper.hpp>
#include <medusa/bits/approximations/Gaussian.hpp>
#include <medusa/bits/approximations/RBFBasis.hpp>
#include <medusa/bits/approximations/Monomials.hpp>
#include <medusa/bits/approximations/ScaleFunction.hpp>
#include <medusa/bits/approximations/WeightFunction.hpp>
#include <medusa/bits/approximations/Operators.hpp>
#include <medusa/bits/types/Vec.hpp>
#include "Eigen/LU"
#include <cmath>

#include "gtest/gtest.h"

namespace mm {

TEST(Approximations, WLSMon1d) {
    WLS<Monomials<Vec1d>, NoWeight<Vec1d>, ScaleToFarthest,
            Eigen::PartialPivLU<Eigen::MatrixXd>> wls(2, {});

    double h = 0.1;
    wls.compute(0.0, {0.0, -h, h});
    Eigen::VectorXd shape = wls.getShape();
    Eigen::VectorXd expected(3); expected << 1.0, 0.0, 0.0;
    ASSERT_EQ(expected.size(), shape.size());
    EXPECT_NEAR(expected[0], shape[0], 1e-15);
    EXPECT_NEAR(expected[1], shape[1], 1e-15);
    EXPECT_NEAR(expected[2], shape[2], 1e-15);

    shape = wls.getShape(Der1<1>(0));
    expected << 0.0, -1.0/2.0/h, 1.0/2.0/h;
    ASSERT_EQ(expected.size(), shape.size());
    EXPECT_NEAR(expected[0], shape[0], 1e-15);
    EXPECT_NEAR(expected[1], shape[1], 1e-15);
    EXPECT_NEAR(expected[2], shape[2], 1e-15);

    shape = wls.getShape(Der2<1>(0));
    expected << -2.0/h/h, 1.0/h/h, 1.0/h/h;
    ASSERT_EQ(expected.size(), shape.size());
    EXPECT_NEAR(expected[0], shape[0], 1e-13);
    EXPECT_NEAR(expected[1], shape[1], 1e-13);
    EXPECT_NEAR(expected[2], shape[2], 1e-13);
}

TEST(Approximations, WLSMon2d) {
    Vec2d point = 0.0;
    double h = 0.123;
    Range<Vec2d> support = {point, {0, h}, {h, 0}, {-h, 0}, {0, -h},
                                   {-h, h}, {h, h}, {-h, -h}, {h, -h}};

    WLS<Monomials<Vec2d>, NoWeight<Vec2d>, ScaleToClosest,
            Eigen::PartialPivLU<Eigen::MatrixXd>> wls(Monomials<Vec2d>::tensorBasis(2), {});

    wls.compute({0.0, 0.0}, support);
    Eigen::VectorXd shape = wls.getShape(Der2<2>(0)) + wls.getShape(Der2<2>(1, 1));
    Eigen::VectorXd expected(9); expected << -4/h/h, 1/h/h, 1/h/h, 1/h/h, 1/h/h, 0, 0, 0, 0;
    EXPECT_EQ(expected, shape);

    shape = wls.getShape(Der1<2>(0));
    expected << 0, 0, 1/(2*h), -1/(2*h), 0, 0, 0, 0, 0;
    EXPECT_EQ(expected, shape);

    shape = wls.getShape(Der1<2>(1));
    expected << 0, 1/(2*h), 0, 0, -1/(2*h), 0, 0, 0, 0;
    EXPECT_EQ(expected, shape);

    shape = wls.getShape(Der2<2>(0, 1));
    expected << 0, 0, 0, 0, 0, -1./4/h/h, 1./4/h/h, 1./4/h/h, -1./4/h/h;
    EXPECT_EQ(expected, shape);
}

TEST(Approximations, WLSMon2dSingular) {
    Vec2d point = 0.0;
    double h = 0.123;
    Range<Vec2d> support = {point, {0, h}, {h, 0}, {-h, 0}, {-h, h}, {h, h},
                            {-h, 2*h}, {0, 2*h}, {2*h, 0}};

    WLS<Monomials<Vec2d>, NoWeight<Vec2d>, ScaleToClosest,
            Eigen::PartialPivLU<Eigen::MatrixXd>> wls(Monomials<Vec2d>::tensorBasis(2), {});

    wls.compute({0.0, 0.0}, support);
    Eigen::VectorXd shape = wls.getShape(Der2<2>(0)) + wls.getShape(Der2<2>(1));
//    EXPECT_FALSE(shape.allFinite());  // different behaviour in new Eigen version?

    shape = wls.getShape(Der1<2>(0));
//    EXPECT_FALSE(shape.allFinite());  // different behaviour in new Eigen version?

    WLS<Monomials<Vec2d>, NoWeight<Vec2d>, ScaleToClosest,
            JacobiSVDWrapper<double>> wls_svd(Monomials<Vec2d>::tensorBasis(2), {});
    wls_svd.compute({0.0, 0.0}, support);

    shape = wls_svd.getShape(Lap<2>());

    Eigen::VectorXd expected(9);
    expected << 1.0/20/h/h, -2/h/h, -1.0/20/h/h, 13.0/20/h/h, 0, 0, 0, 1.0/h/h, 7.0/20/h/h;
    ASSERT_EQ(expected.size(), shape.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_NEAR(expected[i], shape[i], 2e-13);
    }
    shape = wls_svd.getShape(Der1<2>(1));
    expected << -33/(40*h), 2/h, -27/(40*h), -9/(40*h), 0, 0, 0, -1/(2*h), 9/(40*h);
    ASSERT_EQ(expected.size(), shape.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_NEAR(expected[i], shape[i], 2e-14);
    }
}

TEST(Approximations, WLSGauss2d) {
    double h = 0.123;
    double s = 1.5;

    Gaussian<double> g(s);
    WLS<RBFBasis<Gaussian<double>, Vec2d>, NoWeight<Vec2d>, NoScale,
            Eigen::LLT<Eigen::MatrixXd>> wls({9, g}, {});

    Range<Vec2d> support = {{0, 0}, {0, h}, {h, 0}, {0, -h}, {-h, 0},
                            {-h, h}, {-h, -h}, {h, -h}, {h, h}};
    wls.compute({0.0, 0.0}, support);

    double a = -4.*(s*s + h*h/std::pow(std::sinh((h/s)*(h/s)), 2)) / std::pow(s, 4);
    double b = 4.*std::exp(3.*(h/s)*(h/s))*h*h / std::pow(-1+std::exp(2*(h/s)*(h/s)), 2)
               / std::pow(s, 4);

    Eigen::VectorXd shape = wls.getShape(Lap<2>());
    Eigen::VectorXd expected(9); expected << a, b, b, b, b, 0, 0, 0, 0;
    ASSERT_EQ(expected.size(), shape.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_NEAR(expected[i], shape[i], 1e-5);
    }
}

TEST(Approximations, DISABLED_WLSUsageExmaple) {
    /// [WLS usage example]
    WLS<Monomials<Vec2d>, NoWeight<Vec2d>, ScaleToFarthest, Eigen::PartialPivLU<Eigen::MatrixXd>>
        appr(Monomials<Vec2d>::tensorBasis(2));  // full type specification with scale and solver
    std::cout << appr << std::endl;

    WLS<Monomials<Vec2d>> appr2(2);  // using default parameters

    // Local neighbourhood
    double h = 0.1;
    Range<Vec2d> support = {{0, 0}, {0, h}, {h, 0}, {0, -h}, {-h, 0},
                            {-h, h}, {-h, -h}, {h, -h}, {h, h}};

    // compute approximations at point `{0.0, 0.0}`.
    appr.compute({0.0, 0.0}, support);

    // Get info about the computation.
    std::cout << appr.basis() << std::endl;
    std::cout << appr.weight() << std::endl;
    std::cout << appr.center() << std::endl;
    std::cout << appr.scale() << std::endl;
    std::cout << appr.localCoordinates() << std::endl;
    std::cout << appr.getMatrix() << std::endl;
    std::cout << appr.getWeights() << std::endl;
    Eigen::PartialPivLU<Eigen::MatrixXd> solver = appr.solver();

    // Get shape (stencil weights) for approximation of Laplacian.
    Eigen::VectorXd shape = appr.getShape(Lap<2>());

    // Get shape (stencil weights) for approximation of value
    shape = appr.getShape();
    /// [WLS usage example]
    (void) solver;
}

}  // namespace mm
