#include "medusa/bits/approximations/Multiquadric.hpp"

#include "gtest/gtest.h"

namespace mm {

TEST(Approximations, MultiQuadraticEvaluate) {
    double s = 1/0.66;
    Multiquadric<double> mq(s);
    double r1 = 0;
    double r2 = 1.17;
    EXPECT_NEAR(1, mq(r1), 1e-15);
    EXPECT_NEAR(1.2286789653933200, mq(r2), 1e-15);

    EXPECT_NEAR(0.21780000000000000, mq(r1, 1), 1e-15);
    EXPECT_NEAR(0.17726355389365580, mq(r2, 1), 1e-15);

    EXPECT_NEAR(-0.0474368400000000, mq(r1, 2), 1e-15);
    EXPECT_NEAR(-0.025574107170419560, mq(r2, 2), 1e-15);

    EXPECT_NEAR(-0.01047023439346170, mq(r2, 6), 1e-15);
}

TEST(Approximations, DISABLED_MultiQuadraticUsageExample) {
    /// [Multiquadric usage example]
    Multiquadric<double> m(1.324);
    double val = m(4.0);     // value at 4.0 (this is the squared radius)
    double der = m(0.0, 2);  // second derivative at 0.0
    std::cout << m << std::endl;
    /// [Multiquadric usage example]
    std::cout << val << der;  // to not be unused
}

}  // namespace mm
