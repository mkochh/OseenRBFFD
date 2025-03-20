#include "medusa/bits/domains/compute_normal_fwd.hpp"

#include "gtest/gtest.h"

#include <Eigen/Geometry>

namespace mm {

TEST(Domains, NormalTest3D) {
    Eigen::Matrix<double, 3, 2> J;
    J.col(1) = Vec3d(0.5, 0.5, 0);
    J.col(0) = Vec3d(1, 0.5, 0);
    Vec3d normal = -J.col(0).cross(J.col(1));
    normal.normalize();
    Vec3d normal2 = surface_fill_internal::compute_normal(J);
    EXPECT_LT((normal - normal2).norm(), 1e-15);
}

TEST(Domains, NormalTest1Dto3D) {
    Eigen::Matrix<double, 3, 1> J;
    J.col(0) = Vec3d(0.5, 0.5, 0);
    Vec3d normal = surface_fill_internal::compute_normal(J);
    EXPECT_LT(std::abs(normal.dot(J.col(0))), 1e-15);
}

TEST(Domains, NormalTest2D) {
    Eigen::Matrix<double, 2, 1> J;
    J << 0.1, 2.3;
    Vec2d normal(J(1), -J(0)); normal.normalize();
    Vec2d normal2 = surface_fill_internal::compute_normal(J);
    EXPECT_LT((normal - normal2).norm(), 1e-15);
}

}  // namespace mm
