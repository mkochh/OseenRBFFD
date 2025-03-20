#include <medusa/bits/approximations/WeightFunction.hpp>
#include <medusa/bits/approximations/Gaussian.hpp>
#include <medusa/bits/types/Vec_fwd.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(Approximations, RBFWeight2D) {
    double s = 1.23;
    RBFWeight<Gaussian<double>, Vec2d> weight(s);
    Vec2d p(1.37, -0.97);
    EXPECT_NEAR(0.15528149718493840, weight(p), 1e-15);
//    EXPECT_NEAR(-0.28122896575235060, weight(p, {1, 0}), 1e-14);
//    EXPECT_NEAR(0.545584193559560000, weight(p, {1, 1}), 1e-14);
//    EXPECT_NEAR(0.050053899223872840, weight(p, {0, 2}), 1e-12);
}

TEST(Approximations, DISABLED_RBFWeightUsageExample) {
    /// [RBF weight usage example]
    double s = 1.23;
    RBFWeight<Gaussian<double>, Vec2d> weight(s);
    Vec2d p(1.37, -0.97);
    double val = weight(p);
    std::cout << weight << std::endl;
    /// [RBF weight usage example]
    (void) val;
}

}  // namespace mm
