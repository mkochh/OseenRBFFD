#include <medusa/bits/approximations/Gaussian.hpp>
#include <medusa/bits/approximations/RBFFD.hpp>
#include "medusa/bits/approximations/RBFInterpolant.hpp"

#include "gtest/gtest.h"

namespace mm {

TEST(Approximations, RBFInterpolantRBF) {
    double h = 0.154, u1 = 3.23432, u2 = -2.3234, u3 = 0.12443498, u4 = 1.908432, u5 = -0.98742532;
    Vec2d c = {0.05, 0.05};
    std::vector<Vec2d> pts = {{0, 0}, {0, h}, {0, -h}, {h, 0}, {-h, 0}};
    for (auto& x : pts) x += c;

    Eigen::VectorXd values(pts.size()); values << u1, u2, u3, u4, u5;
    RBFFD<Gaussian<double>, Vec2d, ScaleToFarthest> approx(1.0);
    auto appr = approx.getApproximant({0.123, -0.654}, pts, values);
    for (int i = 0; i < 5; ++i) {
        EXPECT_NEAR(values[i], appr(pts[i]), 1e-12);
    }
}

TEST(Approximations, RBFInterpolantRBFPoly) {
    double h = 0.154, u1 = 3.23432, u2 = -2.3234, u3 = 0.12443498, u4 = 1.908432, u5 = -0.98742532;
    Vec2d c = {0.05, 0.05};
    std::vector<Vec2d> pts = {{0, 0}, {0, h}, {0, -h}, {h, 0}, {-h, 0}};
    for (auto& x : pts) x += c;

    Eigen::VectorXd values(pts.size()); values << u1, u2, u3, u4, u5;
    RBFFD<Gaussian<double>, Vec2d, ScaleToFarthest> approx(1.0, 2);
    auto appr = approx.getApproximant({0.123, -0.654}, pts, values);
    for (int i = 0; i < 5; ++i) {
        EXPECT_NEAR(values[i], appr(pts[i]), 1e-13);
    }
    EXPECT_NEAR(appr.coefficients().head(pts.size()).sum(), 0, 1e-13);
}

TEST(Approximations, RBFInterpolantCoeff) {
    std::vector<Vec2d> p = {
            {-0.11766288806517933, -0.9112411814898387}, {0.4962142406313803, 0.9404764517552069},
            {-0.19312682412231208, 0.7528692817356553}, {-0.21512291968760122, -0.5492097324624199},
            {0.44926381204412813, -0.587474226339622}, {-0.025435108077836777, 0.8874183726242024},
            {0.5106615291357914, -0.7659609534980718}};
    Eigen::VectorXd values(p.size());
    values << 0.43897754604943806, 0.598385054145514, 0.9727438713783174, 0.41015894702010325,
              0.07276447441095724, 0.895450228637269, 0.5166314722086436;
    Monomials<Vec2d> mon(1);
    Eigen::VectorXd expected(p.size()+mon.size());  // computed with Mathematica (see rbf_fit.nb)
    expected << -2.5839802672432306, 0.05949468465926759, 0.3070276203705349,
                2.6614384716333603, -9.109317862041848, 0.06397458900968854,
                8.601362763612228, 0.7723570730895989, -0.17108612429083783,
                -0.11758425393124862;

    Vec2d c = {0.1, 0.1};
    RBFFD<Gaussian<double>, Vec2d, ScaleToFarthest> approx(1.0, mon);
    auto appr = approx.getApproximant(c, p, values);
    EXPECT_LT((appr.coefficients() - expected).norm(), 1e-13);
}

TEST(Approximation, RBFInterpolantUsage) {
    /// [RBFInterpolant usage example]
    double h = 0.1;
    std::vector<Vec2d> pts = {{0, 0}, {0, h}, {0, -h}, {h, 0}, {-h, 0}};
    Eigen::VectorXd values(pts.size()); values << 0.1, 0.1, 0.1, 0.2, 0.0;  // linear function in x
    RBFFD<Gaussian<double>, Vec2d, ScaleToFarthest> approx(1.0, 2);
    auto appr = approx.getApproximant({0, 0}, pts, values);
    double val = appr({h/2, 0});
    double lap = appr({h/2, 0}, Lap<2>());
    /// [RBFInterpolant usage example]
    EXPECT_NEAR(0.15, val, 1e-16);
    EXPECT_NEAR(0, lap, 1e-16);
}

}  // namespace mm
