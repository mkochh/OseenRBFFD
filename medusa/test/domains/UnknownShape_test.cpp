#include <medusa/bits/domains/UnknownShape.hpp>
#include <medusa/bits/domains/DomainDiscretization.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(Domains, UnknownShape) {
    UnknownShape<Vec2d> shape;
    EXPECT_TRUE(shape.contains({1.0, -2}));
    EXPECT_TRUE(shape.contains({0.0, 0.0}));
    EXPECT_DEATH(shape.discretizeBoundaryWithStep(0.1), "This function is not available");
    EXPECT_DEATH(shape.bbox(), "This function is not available");
}

}  // namespace mm
