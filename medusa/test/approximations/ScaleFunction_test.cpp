#include <medusa/bits/approximations/ScaleFunction.hpp>
#include <medusa/bits/types/Vec.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(Approximations, ScaleToClosest) {
EXPECT_NEAR(1, ScaleToClosest::scale<Vec2d>({0, 0}, {{0, 0}, {0, 1}, {2, 0}}), 1e-15);
EXPECT_NEAR(1, ScaleToClosest::scale<Vec2d>({0, 0}, {{0, 1}, {2, 0}}), 1e-15);
}

}  // namespace mm
