#include <medusa/bits/domains/ShapeUnion.hpp>
#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/domains/BallShape.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(Domains, ShapeUnionContains1D) {
    BallShape<Vec1d> c(0.0, 1.0);
    BoxShape<Vec1d> b(0.5, 1.5);
    auto u = c + b;

    EXPECT_TRUE(u.contains(1.0));
    EXPECT_TRUE(u.contains(0.0));
    EXPECT_TRUE(u.contains(1.5));
    EXPECT_TRUE(u.contains(-1.0));
    EXPECT_FALSE(u.contains(-3.0));
    EXPECT_FALSE(u.contains(5.0));
}

TEST(Domains, ShapeUnionContains2D) {
    BallShape<Vec2d> c(0.0, 1.0);
    BoxShape<Vec2d> b(0.0, 1.0);
    auto u = c + b;

    EXPECT_TRUE(u.contains(0.0));
    EXPECT_TRUE(u.contains(1.0));
    EXPECT_TRUE(u.contains({-1.0, 0.0}));
    EXPECT_FALSE(u.contains({-3.0, 0.2}));
}

TEST(Domains, ShapeUnionContains3D) {
    BallShape<Vec3d> c(0.0, 1.0);
    BoxShape<Vec3d> b(0.0, 1.0);
    auto u = c + b;

    EXPECT_TRUE(u.contains(0.0));
    EXPECT_TRUE(u.contains(1.0));
    EXPECT_TRUE(u.contains({-1.0, 0.0, 0.0}));
    EXPECT_FALSE(u.contains({-3.0, 0.2, 1.0}));
}

TEST(Domains, ShapeUnionDiscretization) {
    BallShape<Vec2d> c(0.0, 1.0);
    BoxShape<Vec2d> b(0.0, 1.0);
    auto u = c + b;
    (void) dynamic_cast<const BallShape<Vec2d>&>(u.first());
    (void) dynamic_cast<const BoxShape<Vec2d>&>(u.second());

    auto d = u.discretizeWithStep(0.1);
    EXPECT_EQ(356, d.size());
    for (const auto& p : d.positions()) {
        ASSERT_TRUE(d.shape().contains(p));
    }
    // for actual discretization tests see DomainDiscretization::add.
}

TEST(Domains, ShapeUnionUsageExample) {
    /// [ShapeUnion usage example]
    BoxShape<Vec3d> c(-1.0, 1.0);
    BallShape<Vec3d> b(0.0, 1.0);
    auto u = c + b;
    if (u.contains({0, 0.2, -0.3})) {
        // do something
    }
    std::cout << u << std::endl;
    /// [ShapeUnion usage example]
}

}  // namespace mm
