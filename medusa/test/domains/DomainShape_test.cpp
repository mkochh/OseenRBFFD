#include <medusa/bits/domains/DomainShape.hpp>
#include <medusa/bits/domains/BallShape.hpp>
#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/domains/ShapeUnion.hpp>
#include <medusa/bits/domains/ShapeDifference.hpp>
#include "gtest/gtest.h"

namespace mm {

TEST(Domains, ContainsNested2D) {
    // test inclusion - exclusion
    BoxShape<Vec2d> inner({1.5, 1.5}, {2.5, 2.5});
    BallShape<Vec2d> middle({2, 2}, 1);
    BoxShape<Vec2d> outer({0, 0}, {4, 4});
    auto u1 = outer - middle;
    auto u2 = u1 + inner;
    EXPECT_TRUE(u2.contains({1, 1}));
    EXPECT_FALSE(u2.contains({1.2, 2}));
    EXPECT_TRUE(u2.contains({2, 2}));
}

TEST(Domains, ContainsNested3D) {
    // test inclusion - exclusion
    BoxShape<Vec3d> inner({1.5, 1.5, 1.5}, {2.5, 2.5, 2.5});
    BallShape<Vec3d> middle({2, 2, 2}, 1);
    BoxShape<Vec3d> outer({0, 0, 0}, {4, 4, 4});
    auto u1 = outer - middle;
    auto u2 = u1 + inner;
    EXPECT_TRUE(u2.contains({1, 1, 1}));
    EXPECT_FALSE(u2.contains({1.2, 2, 2}));
    EXPECT_TRUE(u2.contains({2, 2, 2}));
}

TEST(Domains, DISABLED_DomainShapeUsageExample) {
    /// [Domain shape usage example]
    BoxShape<Vec3d> inner({1.5, 1.5, 1.5}, {2.5, 2.5, 2.5});
    BallShape<Vec3d> middle({2, 2, 2}, 1);
    BoxShape<Vec3d> outer({0, 0, 0}, {4, 4, 4});
    auto u = outer - middle + inner;
    std::cout << u << std::endl;
    /// [Domain shape usage example]
}

}  // namespace mm
