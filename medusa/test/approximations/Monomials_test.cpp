#include <medusa/bits/approximations/Monomials.hpp>
#include <medusa/bits/approximations/Operators.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(Approximations, MonomialsTensorBaiss) {
    Monomials<Vec2d> basis = Monomials<Vec2d>::tensorBasis(2);

    Vec2d p = {0.3, 2.3};
    // Values
    EXPECT_DOUBLE_EQ(1, basis.eval(0, p));
    EXPECT_DOUBLE_EQ(p[1], basis.eval(1, p));
    EXPECT_DOUBLE_EQ(p[1] * p[1], basis.eval(2, p));
    EXPECT_DOUBLE_EQ(p[0], basis.eval(3, p));
    EXPECT_DOUBLE_EQ(p[0] * p[1], basis.eval(4, p));
    EXPECT_DOUBLE_EQ(p[0] * p[1] * p[1], basis.eval(5, p));
    EXPECT_DOUBLE_EQ(p[0] * p[0], basis.eval(6, p));
    EXPECT_DOUBLE_EQ(p[0] * p[0] * p[1], basis.eval(7, p));
    EXPECT_DOUBLE_EQ(p[0] * p[0] * p[1] * p[1], basis.eval(8, p));
}

TEST(Approximations, Monomials1D) {
    Monomials<Vec1d> basis(9);
    EXPECT_EQ(10, basis.size());

    Range<Vec1d> points = linspace(Vec1d(-1.0), {1.0}, 20);
    for (auto& p : points) {
        // Values
        EXPECT_DOUBLE_EQ(1, basis.eval(0, p));
        EXPECT_DOUBLE_EQ(p[0], basis.eval(1, p));
        EXPECT_DOUBLE_EQ(p[0] * p[0], basis.eval(2, p));
        // First derivatives
        EXPECT_DOUBLE_EQ(0, basis.evalOp(0, p, Der1<1>(0)));
        EXPECT_DOUBLE_EQ(1, basis.evalOp(1, p, Der1<1>(0)));
        EXPECT_DOUBLE_EQ(2 * p[0], basis.evalOp(2, p, Der1<1>(0)));
        // Second derivatives
        EXPECT_DOUBLE_EQ(0, basis.evalOp(0, p, Der2<1>(0)));
        EXPECT_DOUBLE_EQ(0, basis.evalOp(1, p, Der2<1>(0)));
        EXPECT_DOUBLE_EQ(2, basis.evalOp(2, p, Der2<1>(0)));
    }
}

TEST(Approximations, Monomials2DAtZero) {
    Monomials<Vec2d> basis({{2, 3}, {0, 0}, {1, 2}, {1, 0}});
    EXPECT_EQ(0, basis.evalAt0(0));
    EXPECT_EQ(1, basis.evalAt0(1));
    EXPECT_EQ(0, basis.evalAt0(2));
    EXPECT_EQ(0, basis.evalAt0(3));

    EXPECT_EQ(2*3*2, basis.evalOpAt0(0, Derivative<2>({2, 3})));
    EXPECT_EQ(0, basis.evalOpAt0(1, Der1<2>(0)));
    EXPECT_EQ(0, basis.evalOpAt0(2, Derivative<2>({2, 2})));
    EXPECT_EQ(0, basis.evalOpAt0(2, Der2<2>(1)));
    EXPECT_EQ(2, basis.evalOpAt0(2, Derivative<2>({1, 2})));
    EXPECT_EQ(1, basis.evalOpAt0(3, Der1<2>(0)));
    EXPECT_EQ(0, basis.evalOpAt0(3, Der2<2>(0, 1)));
}

template <int dim>
int find_idx(const Eigen::Matrix<int, dim, Eigen::Dynamic>& powers,
             const Eigen::Matrix<int, dim, 1>& mon) {
    int idx = -1;
    for (int i = 0; i < powers.cols(); ++i) {
        EXPECT_FALSE(powers.col(i) == mon && idx != -1);
        if (powers.col(i) == mon) {
            idx = i;
        }
    }
    EXPECT_NE(idx, -1);
    return idx;
}

TEST(Approximations, Monomials2D) {
    Monomials<Vec2d> basis(9);
    EXPECT_EQ(55, basis.size());

    int c = find_idx(basis.powers(), {0, 0});
    int x = find_idx(basis.powers(), {1, 0});
    int y = find_idx(basis.powers(), {0, 1});
    int xx = find_idx(basis.powers(), {2, 0});
    int xy = find_idx(basis.powers(), {1, 1});
    int yy = find_idx(basis.powers(), {0, 2});

    Range<Vec2d> points = linspace(Vec2d(-1.0), Vec2d(1.0), {30, 30});
    for (auto& p : points) {
        // Values
        ASSERT_DOUBLE_EQ(p[1], basis.eval(y, p));
        ASSERT_DOUBLE_EQ(p[0], basis.eval(x, p));
        ASSERT_DOUBLE_EQ(p[1] * p[1], basis.eval(yy, p));
        ASSERT_DOUBLE_EQ(p[0] * p[1], basis.eval(xy, p));
        ASSERT_DOUBLE_EQ(p[0] * p[0], basis.eval(xx, p));
        // x derivatives
        ASSERT_DOUBLE_EQ(0, basis.evalOp(c, p, Der1<2>(0)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(y, p, Der1<2>(0)));
        ASSERT_DOUBLE_EQ(1, basis.evalOp(x, p, Der1<2>(0)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(yy, p, Der1<2>(0)));
        ASSERT_DOUBLE_EQ(p[1], basis.evalOp(xy, p, Der1<2>(0)));
        ASSERT_DOUBLE_EQ(2 * p[0], basis.evalOp(xx, p, Der1<2>(0)));
        // y derivatives
        ASSERT_DOUBLE_EQ(0, basis.evalOp(c, p, Der1<2>(1)));
        ASSERT_DOUBLE_EQ(1, basis.evalOp(y, p, Der1<2>(1)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(x, p, Der1<2>(1)));
        ASSERT_DOUBLE_EQ(2 * p[1], basis.evalOp(yy, p, Der1<2>(1)));
        ASSERT_DOUBLE_EQ(p[0], basis.evalOp(xy, p, Der1<2>(1)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(xx, p, Der1<2>(1)));
        // xx derivatives
        ASSERT_DOUBLE_EQ(0, basis.evalOp(c, p, Der2<2>(0)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(y, p, Der2<2>(0)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(x, p, Der2<2>(0)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(yy, p, Der2<2>(0)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(xy, p, Der2<2>(0)));
        ASSERT_DOUBLE_EQ(2, basis.evalOp(xx, p, Der2<2>(0)));
        // yy derivatives
        ASSERT_DOUBLE_EQ(0, basis.evalOp(c, p, Der2<2>(1)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(y, p, Der2<2>(1)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(x, p, Der2<2>(1)));
        ASSERT_DOUBLE_EQ(2, basis.evalOp(yy, p, Der2<2>(1)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(xy, p, Der2<2>(1)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(xx, p, Der2<2>(1)));
        // xy derivatives
        ASSERT_DOUBLE_EQ(0, basis.evalOp(c, p, Der2<2>(0, 1)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(y, p, Der2<2>(0, 1)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(x, p, Der2<2>(0, 1)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(yy, p, Der2<2>(0, 1)));
        ASSERT_DOUBLE_EQ(1, basis.evalOp(xy, p, Der2<2>(0, 1)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(xx, p, Der2<2>(0, 1)));
        // laplacian
        ASSERT_DOUBLE_EQ(basis.evalOp(c,  p, Der2<2>(0)) + basis.evalOp(c,  p, Der2<2>(1)), basis.evalOp(c,  p, Lap<2>()));  // NOLINT
        ASSERT_DOUBLE_EQ(basis.evalOp(y,  p, Der2<2>(0)) + basis.evalOp(y,  p, Der2<2>(1)), basis.evalOp(y,  p, Lap<2>()));  // NOLINT
        ASSERT_DOUBLE_EQ(basis.evalOp(x,  p, Der2<2>(0)) + basis.evalOp(x,  p, Der2<2>(1)), basis.evalOp(x,  p, Lap<2>()));  // NOLINT
        ASSERT_DOUBLE_EQ(basis.evalOp(yy, p, Der2<2>(0)) + basis.evalOp(yy, p, Der2<2>(1)), basis.evalOp(yy, p, Lap<2>()));  // NOLINT
        ASSERT_DOUBLE_EQ(basis.evalOp(xy, p, Der2<2>(0)) + basis.evalOp(xy, p, Der2<2>(1)), basis.evalOp(xy, p, Lap<2>()));  // NOLINT
        ASSERT_DOUBLE_EQ(basis.evalOp(xx, p, Der2<2>(0)) + basis.evalOp(xx, p, Der2<2>(1)), basis.evalOp(xx, p, Lap<2>()));  // NOLINT
    }
}

TEST(Approximations, Monomials3D) {
    Monomials<Vec3d> basis(4);
    EXPECT_EQ(35, basis.size());

    int c = find_idx(basis.powers(), {0, 0, 0});
    int x = find_idx(basis.powers(), {1, 0, 0});
    int y = find_idx(basis.powers(), {0, 1, 0});
    int z = find_idx(basis.powers(), {0, 0, 1});
    int xx = find_idx(basis.powers(), {2, 0, 0});
    int xy = find_idx(basis.powers(), {1, 1, 0});
    int xz = find_idx(basis.powers(), {1, 0, 1});
    int yy = find_idx(basis.powers(), {0, 2, 0});
    int yz = find_idx(basis.powers(), {0, 1, 1});
    int zz = find_idx(basis.powers(), {0, 0, 2});

    Range<Vec3d> points = linspace(Vec3d(-1.0), Vec3d(1.0), {20, 20, 20});
    for (auto& p : points) {
        // Values
        ASSERT_DOUBLE_EQ(1, basis.eval(c, p));
        ASSERT_DOUBLE_EQ(p[2], basis.eval(z, p));
        ASSERT_DOUBLE_EQ(p[1], basis.eval(y, p));
        ASSERT_DOUBLE_EQ(p[0], basis.eval(x, p));
        ASSERT_DOUBLE_EQ(p[2] * p[2], basis.eval(zz, p));
        ASSERT_DOUBLE_EQ(p[2] * p[1], basis.eval(yz, p));
        ASSERT_DOUBLE_EQ(p[1] * p[1], basis.eval(yy, p));
        ASSERT_DOUBLE_EQ(p[0] * p[2], basis.eval(xz, p));
        ASSERT_DOUBLE_EQ(p[0] * p[1], basis.eval(xy, p));
        ASSERT_DOUBLE_EQ(p[0] * p[0], basis.eval(xx, p));
        // x derivatives
        ASSERT_DOUBLE_EQ(0, basis.evalOp(c, p, Der1<3>(0)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(z, p, Der1<3>(0)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(y, p, Der1<3>(0)));
        ASSERT_DOUBLE_EQ(1, basis.evalOp(x, p, Der1<3>(0)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(zz, p, Der1<3>(0)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(yz, p, Der1<3>(0)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(yy, p, Der1<3>(0)));
        ASSERT_DOUBLE_EQ(p[2], basis.evalOp(xz, p, Der1<3>(0)));
        ASSERT_DOUBLE_EQ(p[1], basis.evalOp(xy, p, Der1<3>(0)));
        ASSERT_DOUBLE_EQ(2 * p[0], basis.evalOp(xx, p, Der1<3>(0)));
        // yz derivatives
        ASSERT_DOUBLE_EQ(0, basis.evalOp(c, p, Der2<3>(1, 2)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(z, p, Der2<3>(1, 2)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(y, p, Der2<3>(1, 2)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(x, p, Der2<3>(1, 2)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(zz, p, Der2<3>(1, 2)));
        ASSERT_DOUBLE_EQ(1, basis.evalOp(yz, p, Der2<3>(1, 2)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(yy, p, Der2<3>(1, 2)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(xz, p, Der2<3>(1, 2)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(xy, p, Der2<3>(1, 2)));
        ASSERT_DOUBLE_EQ(0, basis.evalOp(xx, p, Der2<3>(1, 2)));

        ASSERT_NEAR(basis.evalOp(c,  p, Der2<3>(0)) + basis.evalOp(c,  p, Der2<3>(1)) + basis.evalOp(c,  p, Der2<3>(2)), basis.evalOp(c,  p, Lap<3>()), 1e-14);  // NOLINT
        ASSERT_NEAR(basis.evalOp(z,  p, Der2<3>(0)) + basis.evalOp(z,  p, Der2<3>(1)) + basis.evalOp(z,  p, Der2<3>(2)), basis.evalOp(z,  p, Lap<3>()), 1e-14);  // NOLINT
        ASSERT_NEAR(basis.evalOp(y,  p, Der2<3>(0)) + basis.evalOp(y,  p, Der2<3>(1)) + basis.evalOp(y,  p, Der2<3>(2)), basis.evalOp(y,  p, Lap<3>()), 1e-14);  // NOLINT
        ASSERT_NEAR(basis.evalOp(x,  p, Der2<3>(0)) + basis.evalOp(x,  p, Der2<3>(1)) + basis.evalOp(x,  p, Der2<3>(2)), basis.evalOp(x,  p, Lap<3>()), 1e-14);  // NOLINT
        ASSERT_NEAR(basis.evalOp(zz, p, Der2<3>(0)) + basis.evalOp(zz, p, Der2<3>(1)) + basis.evalOp(zz, p, Der2<3>(2)), basis.evalOp(zz, p, Lap<3>()), 1e-14);  // NOLINT
        ASSERT_NEAR(basis.evalOp(yz, p, Der2<3>(0)) + basis.evalOp(yz, p, Der2<3>(1)) + basis.evalOp(yz, p, Der2<3>(2)), basis.evalOp(yz, p, Lap<3>()), 1e-14);  // NOLINT
        ASSERT_NEAR(basis.evalOp(yy, p, Der2<3>(0)) + basis.evalOp(yy, p, Der2<3>(1)) + basis.evalOp(yy, p, Der2<3>(2)), basis.evalOp(yy, p, Lap<3>()), 1e-14);  // NOLINT
        ASSERT_NEAR(basis.evalOp(xz, p, Der2<3>(0)) + basis.evalOp(xz, p, Der2<3>(1)) + basis.evalOp(xz, p, Der2<3>(2)), basis.evalOp(xz, p, Lap<3>()), 1e-14);  // NOLINT
        ASSERT_NEAR(basis.evalOp(xy, p, Der2<3>(0)) + basis.evalOp(xy, p, Der2<3>(1)) + basis.evalOp(xy, p, Der2<3>(2)), basis.evalOp(xy, p, Lap<3>()), 1e-14);  // NOLINT
        ASSERT_NEAR(basis.evalOp(xx, p, Der2<3>(0)) + basis.evalOp(xx, p, Der2<3>(1)) + basis.evalOp(xx, p, Der2<3>(2)), basis.evalOp(xx, p, Lap<3>()), 1e-14);  // NOLINT
    }
}

TEST(Approximations, MonomialsManual) {
    Monomials<Vec2d> basis({{0, 2}, {1, 2}, {3, 1}, {0, 0}});
    auto points = linspace(Vec2d(-1.0), Vec2d(1.0), {20, 20});
    for (const auto& p : points) {
        EXPECT_DOUBLE_EQ(p[1] * p[1], basis.eval(0, p));
        EXPECT_DOUBLE_EQ(p[0] * p[1] * p[1], basis.eval(1, p));
        EXPECT_DOUBLE_EQ(p[0] * p[0] * p[0] * p[1], basis.eval(2, p));
        EXPECT_DOUBLE_EQ(1, basis.eval(3, p));
    }
}

TEST(Approximations, MonomialsEmpty) {
    Monomials<Vec2d> empty;
    EXPECT_EQ(0, empty.size());
    Monomials<Vec1d> empty1(-1);
    EXPECT_EQ(0, empty1.size());
    EXPECT_EQ(0, Monomials<Vec3d>::tensorBasis(-1).size());
}

TEST(Appoximations, DISABLED_MonomialsUsageExample) {
    /// [Monomials usage example]
    Monomials<Vec2d> basis(2);  // {1, y, x, y^2, xy, x^2}
    std::cout << basis << std::endl;  // Monomials 2D: [0, 0; 0, 1; 1, 0; 0, 2; 1, 1; 2, 0], ...
    double value = basis.eval(2, {0.3, -2.4});  // value = 0.3

    basis = Monomials<Vec2d>::tensorBasis(1);  // tensor basis {1, x, y, xy}
    basis = Monomials<Vec2d>({{1, 0}, {2, 3}});  // from powers
    // Evaluate derivative or any other supported operator.
    double der = basis.evalOp(1, {0.3, -2.4}, Derivative<2>({2, 1}));  // der = 2*3*(-2.4)^2
    /// [Monomials usage example]
    EXPECT_EQ(2*3*(-2.4)*(-2.4), der);
    EXPECT_EQ(0.3, value);
}

}  // namespace mm
