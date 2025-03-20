#include <medusa/bits/domains/ShapeDifference.hpp>

#include "gtest/gtest.h"
#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/domains/BallShape.hpp>
#include <medusa/bits/io/HDF.hpp>

namespace mm {

TEST(Domains, ShapeDifferenceContains1D) {
    BallShape<Vec1d> c(0.0, 1.5);
    BoxShape<Vec1d> b(0.0, 1.5);
    auto u = c - b;

    EXPECT_TRUE(u.contains(0.0));
    EXPECT_TRUE(u.contains(-0.5));
    EXPECT_TRUE(u.contains(-1.5));
    EXPECT_TRUE(u.contains(1.5));  // boundary points are not removed
    EXPECT_FALSE(u.contains(-2.0));
}

TEST(Domains, ShapeDifferenceContains2D) {
    BoxShape<Vec2d> c(-1.0, 1.0);
    BallShape<Vec2d> b(0.0, 1.0);
    auto u = c - b;

    EXPECT_TRUE(u.contains(-1.0));
    EXPECT_TRUE(u.contains(1.0));
    EXPECT_TRUE(u.contains({-1.0, 0.0}));
    EXPECT_FALSE(u.contains(0.0));
    EXPECT_FALSE(u.contains({0.0, 0.9}));
}

TEST(Domains, ShapeDifferenceContains3D) {
    BoxShape<Vec3d> c(-1.0, 1.0);
    BallShape<Vec3d> b(0.0, 1.0);
    auto u = c - b;

    EXPECT_TRUE(u.contains(-1.0));
    EXPECT_TRUE(u.contains(1.0));
    EXPECT_TRUE(u.contains({-1.0, 0.0, 0.0}));
    EXPECT_FALSE(u.contains(0.0));
    EXPECT_FALSE(u.contains({0.0, 0.9, 0.0}));
}

TEST(Domains, ShapeDifferenceDiscretization) {
    BallShape<Vec2d> c(0.0, 1.0);
    BoxShape<Vec2d> b(0.0, 1.1);
    auto u = c - b;
    (void) dynamic_cast<const BallShape<Vec2d>&>(u.first());
    (void) dynamic_cast<const BoxShape<Vec2d>&>(u.second());

    auto d = u.discretizeWithStep(0.1);
    EXPECT_EQ(269, d.size());
    for (const auto& p : d.positions()) {
        ASSERT_TRUE(d.shape().contains(p)) << p;
    }
    // for actual discretization tests see DomainDiscretization::subtract.
}

TEST(Domains, ShapeDifferenceUsageExample) {
    /// [ShapeDifference usage example]
    BoxShape<Vec3d> c(-1.0, 1.0);
    BallShape<Vec3d> b(0.0, 1.0);
    auto u = c - b;
    if (u.contains({0, 0.2, -0.3})) {
        // do something
    }
    std::cout << u << std::endl;
    /// [ShapeDifference usage example]
}

}  // namespace mm
