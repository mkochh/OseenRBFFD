#include "medusa/bits/approximations/Polyharmonic.hpp"
#include <cmath>

#include "gtest/gtest.h"

namespace mm {
/*
TEST(Approximations, PolyharmonicEvaluate) {
    Polyharmonic<double> ph(1);
    double r0 = 0;
    double r2 = 1.23;
    EXPECT_NEAR(0, ph(r0), 1e-15);
    EXPECT_NEAR(std::sqrt(r2), ph(r2), 1e-15);

    EXPECT_NEAR(0, ph(r0, 1), 1e-15);
    EXPECT_NEAR(0.5 / std::sqrt(r2), ph(r2, 1), 1e-15);

    EXPECT_NEAR(0, ph(r0, 2), 1e-15);
    EXPECT_NEAR(-0.25 / r2 / std::sqrt(r2), ph(r2, 2), 1e-15);
}
*/

TEST(Approximations, PolyharmonicEvaluateCompileTime) {
    Polyharmonic<double, 3> ph;
    double r0 = 0;
    double r2 = 1.23;
    EXPECT_NEAR(0, ph(r0), 1e-15);
    EXPECT_NEAR(std::sqrt(r2*r2*r2), ph(r2), 1e-15);

    EXPECT_NEAR(0, ph(r0, 0), 1e-15);
    EXPECT_NEAR(std::sqrt(r2*r2*r2), ph(r2), 1e-15);

    EXPECT_NEAR(0, ph(r0, 1), 1e-15);
    EXPECT_NEAR(1.5*std::sqrt(r2), ph(r2, 1), 1e-15);

    EXPECT_NEAR(0, ph(r0, 2), 1e-15);
    EXPECT_NEAR(0.75/std::sqrt(r2), ph(r2, 2), 1e-15);

    EXPECT_NEAR(0, ph(r0, 3), 1e-15);
    EXPECT_NEAR(-0.375 / r2 / std::sqrt(r2), ph(r2, 3), 1e-15);
}

TEST(Approximations, DISABLED_PolyharmonicUsageExample) {
    /// [Polyharmonic usage example]
    Polyharmonic<double, 3> ph;  // order given compile time
    double val = ph(4.0);     // value at 4.0 (this is the squared radius)
    double der = ph(0.0, 2);  // second derivative at 0.0
    std::cout << ph << std::endl;
    /// [Polyharmonic usage example]
    std::cout << val << der;  // to not be unused
}

}  // namespace mm
