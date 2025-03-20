#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/domains/BallShape.hpp>
#include "medusa/bits/domains/RotatedShape.hpp"
#include <Eigen/Geometry>

#include "gtest/gtest.h"

namespace mm {

TEST(Domains, RotatedContains) {
    BoxShape<Vec2d> box2d({-1.1, -1.1}, {1.1, 1.1});
    Eigen::Rotation2Dd Q(PI / 4);
    auto t = box2d.rotate(Q.toRotationMatrix());

    EXPECT_TRUE(t.contains({0, 0}));
    EXPECT_FALSE(t.contains({-1, -1}));
    EXPECT_FALSE(t.contains({1, -1}));
    EXPECT_FALSE(t.contains({-1, 1}));
    EXPECT_FALSE(t.contains({1, 1}));
    EXPECT_TRUE(t.contains({0, 1.3}));
    EXPECT_TRUE(t.contains({0, -1.3}));
    EXPECT_TRUE(t.contains({-1.3, 0}));
    EXPECT_TRUE(t.contains({1.3, 0}));
}

TEST(Domains, RotatedBbox) {
    BoxShape<Vec2d> box2d({-1.1, -1.1}, {1.1, 1.1});
    Eigen::Rotation2Dd Q(PI/4);
    auto t = box2d.rotate(Q.toRotationMatrix());
    double bbox_width = 1.1 * std::sqrt(2);
    auto bbox = t.bbox();
    EXPECT_LT((Vec2d(-bbox_width, -bbox_width) - bbox.first).norm(), 1e-15);
    EXPECT_LT((Vec2d(bbox_width, bbox_width) - bbox.second).norm(), 1e-15);
}


TEST(Domains, RotatedPointsAndNormals) {
    BoxShape<Vec2d> box2d({-1.1, -1.1}, {1.1, 1.1});
    Eigen::Rotation2Dd Q(PI/4);
    auto d = box2d.discretizeWithStep(0.1);
    auto d_rotated = box2d.discretizeWithStep(0.1).rotate(Q.toRotationMatrix());

    for (int i = 0; i < d.size(); ++i) {
        EXPECT_EQ(Q*d.pos(i), d_rotated.pos(i));
        if (d.type(i) < 0) {
            EXPECT_EQ(Q * d.normal(i), d_rotated.normal(i));
        }
    }

    // test 2D overloads
    Eigen::Matrix2d M = Q.toRotationMatrix();
    d.rotate(PI/4);
    auto sh = dynamic_cast<const RotatedShape<Vec2d>*>(&d.shape());
    Eigen::Matrix2d Qd = sh->rotation();
    Eigen::Matrix2d Qd2 = box2d.rotate(PI/4).rotation();
    EXPECT_LT((M-Qd).norm(), 1e-15);
    EXPECT_LT((M-Qd2).norm(), 1e-15);
}

TEST(Domains, RotatedShapeCollapse) {
    BoxShape<Vec3d> box({0, 0, 0}, {1, 1, 1});
    Eigen::Matrix3d Q1 = Eigen::AngleAxisd(PI/2, Vec3d::UnitZ()).toRotationMatrix();
    Eigen::Matrix3d Q2 = Eigen::AngleAxisd(PI/2, Vec3d::UnitX()).toRotationMatrix();
    auto t = box.rotate(Q1).rotate(Q2);
    EXPECT_NE(nullptr, dynamic_cast<const BoxShape<Vec3d>*>(&t.shape()));
    EXPECT_LT((Q2*Q1-t.rotation()).norm(), 1e-15);
}

TEST(Domains, RotatedDiscretize) {
    /// [RotatedShape usage example]
    Eigen::AngleAxisd Q(PI/6, Vec3d(1, 1, 1).normalized());

    BallShape<Vec3d> ball(0, 1);
    auto d = ball.discretizeBoundaryWithStep(0.1).rotate(Q.toRotationMatrix());
    auto ball2 = ball.rotate(Q.toRotationMatrix());
    auto d2 = ball2.discretizeBoundaryWithStep(0.1);
    std::cout << ball2 << std::endl;
    /// [RotatedShape usage example]

    EXPECT_EQ(d2.positions(), d.positions());
    EXPECT_EQ(d2.types(), d.types());
    EXPECT_EQ(d2.normals(), d.normals());

    d.assert_is_valid();
    d2.assert_is_valid();
}

}  // namespace mm
