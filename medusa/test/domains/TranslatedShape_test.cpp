#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/domains/BallShape.hpp>
#include "medusa/bits/domains/TranslatedShape.hpp"

#include "gtest/gtest.h"

namespace mm {

TEST(Domains, TranslatedContains) {
    BoxShape<Vec2d> box2d({1, 3}, {1.2, 3.5});
    auto t = box2d.translate({-1.1, -3.1});

    EXPECT_TRUE(t.contains({0, 0}));
    EXPECT_TRUE(t.contains({0.1, 0.4}));
    EXPECT_FALSE(t.contains({1.1, 3.1}));
}

TEST(Domains, TranslatedBbox) {
    BoxShape<Vec2d> box2d({-1.1, -1.1}, {1.1, 1.1});
    Vec2d a = {-1, -2};
    auto t = box2d.translate(a);
    auto bbox = t.bbox();
    EXPECT_LT((box2d.beg()+a - bbox.first).norm(), 1e-15);
    EXPECT_LT((box2d.end()+a - bbox.second).norm(), 1e-15);
}

TEST(Domains, TranslatedShapeCollapse) {
    BoxShape<Vec2d> box2d({1, 3}, {1.2, 3.5});
    auto t = box2d.translate({1, 1}).translate({2, 2});
    EXPECT_NE(nullptr, dynamic_cast<const BoxShape<Vec2d>*>(&t.shape()));
    EXPECT_EQ(Vec2d(3, 3), t.translation());
}

TEST(Domains, TranslatedDiscretize) {
    /// [TranslatedShape usage example]
    BallShape<Vec3d> ball(0, 1);
    auto d = ball.discretizeBoundaryWithStep(0.1).translate(1.0);

    auto ball2 = ball.translate(1.0);
    auto d2 = ball2.discretizeBoundaryWithStep(0.1);
    std::cout << ball2 << std::endl;
    /// [TranslatedShape usage example]

    EXPECT_EQ(d2.positions(), d.positions());
    EXPECT_EQ(d2.types(), d.types());
    EXPECT_EQ(d2.normals(), d.normals());

    d.assert_is_valid();
    d2.assert_is_valid();
}

}  // namespace mm
