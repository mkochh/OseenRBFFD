#include <medusa/bits/approximations/Gaussian.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(Approximations, GaussiansEvaluate) {
    double s = 1.43;
    Gaussian<double> g(s);
    double r0 = 0;
    double r2 = 1.23;
    EXPECT_NEAR(1, g(r0), 1e-15);
    EXPECT_NEAR(0.54799100538918500839, g(r2), 1e-15);

    EXPECT_NEAR(-1/s/s, g(r0, 1), 1e-15);
    EXPECT_NEAR(-0.26797936592947577309, g(r2, 1), 1e-15);

    EXPECT_NEAR(1/s/s/s/s, g(r0, 2), 1e-15);
    EXPECT_NEAR(0.13104766293191636417, g(r2, 2), 1e-15);
}

TEST(Approximations, DISABLED_GaussianUsageExample) {
    /// [Gaussian usage example]
    Gaussian<double> g(1.324);
    double val = g(4.0);     // value at 4.0 (this is the squared radius)
    double der = g(0.0, 2);  // second derivative at 0.0
    std::cout << g << std::endl;
    /// [Gaussian usage example]
    std::cout << val << der;  // to not be unused
}

}  // namespace mm
