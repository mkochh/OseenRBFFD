#include <medusa/bits/types/Vec.hpp>

#include "gtest/gtest.h"

#include <vector>
#include <algorithm>
#include <medusa/bits/utils/assert.hpp>

namespace mm {

TEST(Types, Matrix) {
    // test for different constructors of `int` fixed size vectors, where constructor call might
    // be ambiguous.
    Eigen::Matrix2d M(2.4);
    EXPECT_EQ(2.4, M(0, 0));
    EXPECT_EQ(2.4, M(0, 1));
    EXPECT_EQ(2.4, M(1, 1));
    EXPECT_EQ(2.4, M(1, 0));

    Eigen::VectorXi v(3);
    v[0] = v[1] = v[2] = -5;
    EXPECT_EQ(-5, v[0]);
    EXPECT_EQ(-5, v[1]);
    EXPECT_EQ(-5, v[2]);

    Eigen::Vector3i v2(6);
    EXPECT_EQ(6, v2[0]);
    EXPECT_EQ(6, v2[1]);
    EXPECT_EQ(6, v2[2]);

    Eigen::VectorXd v3(3);
    EXPECT_EQ(v3.size(), 3);
}

TEST(Types, VecConstruct) {
    Vec1d a(6.0);
    EXPECT_EQ(6.0, a[0]);

    Vec2d b({2.0, 3.0});
    EXPECT_EQ(2.0, b[0]);
    EXPECT_EQ(3.0, b[1]);

    Vec3d c(5.0);
    EXPECT_EQ(5.0, c[0]);
    EXPECT_EQ(5.0, c[1]);
    EXPECT_EQ(5.0, c[2]);

    Vec3d d;
    d << 1, 5, 7;
    EXPECT_EQ(1, d[0]);
    EXPECT_EQ(5, d[1]);
    EXPECT_EQ(7, d[2]);

    Vec2d e(3.4, -1.2);
    EXPECT_EQ(3.4, e[0]);
    EXPECT_EQ(-1.2, e[1]);
}

TEST(Types, VecAssign) {
    Vec2d a;
    a = {1, 2};
    EXPECT_DOUBLE_EQ(1.0, a[0]);
    EXPECT_DOUBLE_EQ(2.0, a[1]);

    a = 5.3;
    EXPECT_EQ(5.3, a[0]);
    EXPECT_EQ(5.3, a[1]);

    Vec2d b, c;
    b = a;
    EXPECT_EQ(a[0], b[0]);
    EXPECT_EQ(a[1], b[1]);

    c = std::move(a);
    EXPECT_EQ(b[0], c[0]);
    EXPECT_EQ(b[1], c[1]);

    a[1] = 4;
    EXPECT_EQ(a[1], 4);
}

TEST(Types, VecCompare) {
    Vec3d a(0.0);
    Vec3d b(0.0);
    EXPECT_TRUE(a == b);
    EXPECT_TRUE(a >= b);
    EXPECT_TRUE(a <= b);
    EXPECT_FALSE(a != b);
    EXPECT_FALSE(a < b);
    EXPECT_FALSE(a > b);
    b[0] = 1;
    EXPECT_FALSE(a == b);
    EXPECT_FALSE(a >= b);
    EXPECT_TRUE(a <= b);
    EXPECT_TRUE(a != b);
    EXPECT_TRUE(a < b);
    EXPECT_FALSE(a > b);
    a = {1, 1, 0};
    b = {1, 1, 1};
    EXPECT_FALSE(a == b);
    EXPECT_FALSE(a >= b);
    EXPECT_TRUE(a <= b);
    EXPECT_TRUE(a != b);
    EXPECT_TRUE(a < b);
    EXPECT_FALSE(a > b);

    EXPECT_LT(Vec2d({0, 0.5}), Vec2d({0.5, 0}));
    EXPECT_FALSE(Vec2d({0.5, 0}) < Vec2d({0, 0.5}));
    EXPECT_GT(Vec2d({0.5, 0}), Vec2d({0, 0.5}));
    EXPECT_FALSE(Vec2d({0, 0.5}) > Vec2d({0.5, 0}));

    std::vector<Vec2d> expected = {{0, 0},   {0, 0.5}, {0, 1},   {0.5, 0},
                                   {0.5, 1}, {1, 0},   {1, 0.5}, {1, 1}};

    for (size_t i = 0; i < expected.size(); ++i) {
        for (size_t j = i + 1; j < expected.size(); ++j) {
            EXPECT_LT(expected[i], expected[j])
                << expected[i].transpose() << " < " << expected[j].transpose();
            EXPECT_GT(expected[j], expected[i])
                << expected[i].transpose() << " > " << expected[j].transpose();
        }
    }

    std::vector<Vec2d> tosort = expected;
    std::sort(tosort.begin(), tosort.end());
    ASSERT_EQ(expected.size(), tosort.size());
    for (size_t i = 0; i < expected.size(); ++i)
        EXPECT_EQ(expected[i], tosort[i])
            << expected[i].transpose() << " == " << tosort[i].transpose();
}

TEST(Types, VecBasicMath) {
    Vec3d a({1, 3, -1.2}), b({1, 0, 2.5});

    EXPECT_EQ(Vec3d({-1, -3, 1.2}), -a);
    EXPECT_EQ(Vec3d({2, 3, 1.3}), a + b);
    EXPECT_EQ(Vec3d({0, 3, -3.7}), a - b);
    EXPECT_EQ(Vec3d({0, -3, 3.7}), b - a);
    Vec3d r = 2.1 * a;
    EXPECT_DOUBLE_EQ(2.1, r[0]);
    EXPECT_DOUBLE_EQ(6.3, r[1]);
    EXPECT_DOUBLE_EQ(-1.2 * 2.1, r[2]);
    r = a * 2.1;
    EXPECT_DOUBLE_EQ(2.1, r[0]);
    EXPECT_DOUBLE_EQ(6.3, r[1]);
    EXPECT_DOUBLE_EQ(-1.2 * 2.1, r[2]);
    r = r / 2.1;
    EXPECT_DOUBLE_EQ(1, r[0]);
    EXPECT_DOUBLE_EQ(3, r[1]);
    EXPECT_DOUBLE_EQ(-1.2, r[2]);
    EXPECT_DOUBLE_EQ(7.25, b.squaredNorm());
    EXPECT_DOUBLE_EQ(2.6925824035672519, b.norm());
    a += b;
    ASSERT_EQ((Vec3d{2, 3, 1.3}), a);
    a *= 2.1;
    EXPECT_DOUBLE_EQ(4.2, a[0]);
    EXPECT_DOUBLE_EQ(6.3, a[1]);
    EXPECT_DOUBLE_EQ(1.3 * 2.1, a[2]);
    a /= 0.5;
    EXPECT_DOUBLE_EQ(8.4, a[0]);
    EXPECT_DOUBLE_EQ(12.6, a[1]);
    EXPECT_DOUBLE_EQ(2 * 1.3 * 2.1, a[2]);
}

TEST(Types, VecIterate) {
    Vec3d a(0.0);
    for (auto& x : a) x += 4;
    EXPECT_EQ((Vec3d{4, 4, 4}), a);
}

TEST(Types, DISABLED_VecUsageExample) {
    /// [Vec usage example]
    Vec3d a(4.0);
    a = {3, -4.3, 2};
    a = 12.3;
    a /= 23.4;
    Vec3d b(12.3, 34.5, -2.1);
    b -= a;
    std::cout << a << std::endl;
    /// [Vec usage example]
}

}  // namespace mm
