#include "medusa/bits/approximations/WLSApproximant.hpp"
#include "medusa/bits/approximations/Operators.hpp"
#include "medusa/bits/approximations/WLS.hpp"
#include "medusa/bits/approximations/JacobiSVDWrapper.hpp"
#include "medusa/bits/approximations/Monomials.hpp"
#include "medusa/bits/approximations/RBFBasis.hpp"
#include "medusa/bits/approximations/WeightFunction.hpp"
#include "medusa/bits/approximations/ScaleFunction.hpp"
#include "medusa/bits/approximations/Gaussian.hpp"

#include "medusa/bits/types/Vec.hpp"

#include "gtest/gtest.h"

namespace mm {

TEST(Approximations, WLSApproximantInterpolation) {
    Vec1d point = 1.2;
    std::vector<Vec1d> support = {0.5, 0.6, 0.8, 1.0, 1.1, 1.4, 1.5};
    int n = support.size();
    Eigen::VectorXd values(n);  // monomial 2x^2 - 3x + 5 is reproduced exactly
    for (int i = 0; i < n; ++i) {
        values[i] = 2*support[i][0]*support[i][0] - 3*support[i][0] + 4;
    }

    WLS<Monomials<Vec1d>, NoWeight<Vec1d>, ScaleToClosest> wls(n-1);
    auto appr = wls.getApproximant(point, support, values);
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(2*support[i][0]*support[i][0] - 3*support[i][0] + 4, appr(support[i]), 1e-13);
        EXPECT_NEAR(4*support[i][0] - 3, appr(support[i], Der1<1>(0)), 1e-13);
        EXPECT_NEAR(4, appr(support[i], Der2<1>(0)), 1e-11);
    }
    EXPECT_NEAR(appr.residual(), 0, 1e-16);
}

TEST(Approximations, WLSApproximantLS) {
    double h = 0.1, u1 = 3.23432, u2 = -2.3234, u3 = 0.12443498, u4 = 1.908432, u5 = -0.98742532;
    std::vector<Vec2d> pts = {{0, 0}, {0, -h}, {0, h}, {h, 0}, {-h, 0}};
    Eigen::VectorXd values(pts.size()); values << u1, u2, u3, u4, u5;
    WLS<Monomials<Vec2d>, NoWeight<Vec2d>, ScaleToFarthest> wls(1);

    Eigen::VectorXd expected_coeff(wls.basis().size());
    expected_coeff << u1/5. + u2/5. + u3/5. + u4/5. + u5/5.,
            -u2/2. + u3/2.,
            u4/2. - u5/2.;

    auto appr = wls.getApproximant(0.0, pts, values);
    Eigen::VectorXd coeff = appr.coefficients();

    for (int i = 0; i < appr.basis().size(); ++i) {
        EXPECT_NEAR(coeff[i], expected_coeff[i], 1e-15);
    }
}

TEST(Approximations, WLSApproximantWLS) {
    /// [WLSApproximant usage example]
    double h = 0.1, u1 = 3.23432, u2 = -2.3234, u3 = 0.12443498, u4 = 1.908432, u5 = -0.98742532;
    Vec2d c = {0.05, 0.05};
    std::vector<Vec2d> pts = {{0, 0}, {0, h}, {0, -h}, {h, 0}, {-h, 0}};
    Eigen::VectorXd values(pts.size()); values << u1, u2, u3, u4, u5;
    WLS<Monomials<Vec2d>, GaussianWeight<Vec2d>, ScaleToClosest> wls(1, 1.0);
    auto appr = wls.getApproximant(c, pts, values);
    double value = appr(0);
    double der = appr(0, Der1<2>(1));
    /// [WLSApproximant usage example]

    EXPECT_NEAR(value, 3.2248113457636491, 1e-15);
    EXPECT_NEAR(der, -55.453110407462376, 1e-15);

    Eigen::VectorXd expected_coeff(wls.basis().size());
    expected_coeff <<  // calculated using Mathematica (see test/approximations/wls_fit.nb)
            0.000668905719969*u1 + 0.499832773570008*u2 - 0.0001672264299923161*u3 +
            0.499832773570008*u4 - 0.0001672264299923161*u5,
            -0.7047418423337238*u1 + 0.7063969822976793*u2 - 0.000709798888868241*u3 -
            0.0004726705375436*u4 - 0.0004726705375436009*u5,
            -0.7047418423337238*u1 - 0.0004726705375436*u2 - 0.0004726705375436009*u3 +
            0.7063969822976793*u4 - 0.000709798888868241*u5;

    Eigen::VectorXd coeff = appr.coefficients();
    for (int i = 0; i < appr.basis().size(); ++i) {
        EXPECT_NEAR(coeff[i], expected_coeff[i], 1e-14);
    }
    EXPECT_NEAR(appr(0), 3.2248113457636491, 1e-15);
}

TEST(Approximations, WLSInterpolantRBF) {
    double h = 0.154, u1 = 3.23432, u2 = -2.3234, u3 = 0.12443498, u4 = 1.908432, u5 = -0.98742532;
    Vec2d c = {0.05, 0.05};
    std::vector<Vec2d> pts = {{0, 0}, {0, h}, {0, -h}, {h, 0}, {-h, 0}};
    for (auto& x : pts) x += c;

    Eigen::VectorXd values(pts.size()); values << u1, u2, u3, u4, u5;
    WLS<Gaussians<Vec2d>, NoWeight<Vec2d>, ScaleToFarthest> wls({5, 1.0});
    auto appr = wls.getApproximant({0.123, -0.654}, pts, values);
    for (int i = 0; i < 5; ++i) {
        EXPECT_NEAR(values[i], appr(pts[i]), 1e-11);
    }
}

}  // namespace mm
