#include "medusa/bits/approximations/InverseMultiquadric.hpp"

#include "gtest/gtest.h"

namespace mm {

TEST(Approximations, InverseMultiQuadraticEvaluate) {
    double s = 1/1.43;
    InverseMultiquadric<double> imq(s);
    double r0 = 0.0;
    double r2 = 1.74;
    EXPECT_NEAR(1, imq(r0), 1e-15);
    EXPECT_NEAR(0.4683891614366080000, imq(r2), 1e-15);

    EXPECT_NEAR(-1.022450000000000000, imq(r0, 1), 1e-15);
    EXPECT_NEAR(-0.105066094730786300, imq(r2, 1), 1e-15);

    EXPECT_NEAR(3.1362120074999990000, imq(r0, 2), 1e-15);
    EXPECT_NEAR(0.0707032858838209400, imq(r2, 2), 1e-15);

    EXPECT_NEAR(0.6202486805208258000, imq(r2, 6), 1e-15);
}

TEST(Approximations, DISABLED_InverseMultiQuadraticUsageExample) {
    /// [InverseMultiquadric usage example]
    InverseMultiquadric<double> im(1.324);
    double val = im(4.0);     // value at 4.0 (this is the squared radius)
    double der = im(0.0, 2);  // second derivative at 0.0
    std::cout << im << std::endl;
    /// [InverseMultiquadric usage example]
    std::cout << val << der;  // to not be unused
}

}  // namespace mm
