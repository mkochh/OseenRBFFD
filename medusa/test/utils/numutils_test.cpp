#include <medusa/bits/utils/numutils.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(Utils, NumSignumZero) {
    assert(signum(0.0f) == 0);
    assert(signum(0l) == 0);
    assert(signum(0ull) == 0);
    assert(signum(0.0l) == 0);
    assert(signum(0) == 0);
    assert(signum(0.0) == 0);
}

TEST(Utils, NumSignumInfNan) {
    // just to define behaviour
    double nan = std::numeric_limits<double>::quiet_NaN();
    double inf = std::numeric_limits<double>::infinity();
    assert(signum(nan) == 0);
    assert(signum(inf) == 1);
    assert(signum(-inf) == -1);
}

TEST(Utils, NumSignumPositive) {
    assert(signum('a') == 1);
    assert(signum(5.4) == 1);
    assert(signum(1e-7) == 1);
    assert(signum(23452345ull) == 1);
    assert(signum(1e20f) == 1);
    assert(signum(-23452345ull) == 1);
}

TEST(Utils, NumSignumNegative) {
    assert(signum(-5.4) == -1);
    assert(signum(-1e-7) == -1);
    assert(signum(-1e20f) == -1);
    assert(signum(-23452345ll) == -1);
}

TEST(Utils, NumFloorCeil) {
    EXPECT_EQ(1, ifloor(1.3));
    EXPECT_EQ(-2, ifloor(-1.3));
    EXPECT_EQ(2, iceil(1.3));
    EXPECT_EQ(-1, iceil(-1.3));
}

TEST(Utils, NumIPowCompileTime) {
    EXPECT_EQ(1, ipow<0>(3.4));
    EXPECT_EQ(1, ipow<0>(-1.4));
    EXPECT_EQ(3.4, ipow<1>(3.4));
    EXPECT_EQ(-1.4, ipow<1>(-1.4));
    EXPECT_EQ(3.4*3.4*3.4, ipow<3>(3.4));
    EXPECT_EQ(1.4*1.4, ipow<2>(-1.4));
}

TEST(Utils, NumIPow) {
    EXPECT_EQ(1, ipow(3.4, 0));
    EXPECT_EQ(1, ipow(-1.4, 0));
    EXPECT_EQ(3.4, ipow(3.4, 1));
    EXPECT_EQ(-1.4, ipow(-1.4, 1));
    EXPECT_EQ(3.4*3.4*3.4, ipow(3.4, 3));
    EXPECT_EQ(1.4*1.4, ipow(-1.4, 2));

    float x = 1.4f;
    float r = 1.0f;
    for (int i = 0; i < 10; ++i) {
        EXPECT_FLOAT_EQ(r, ipow(x, i));
        r *= x;
    }
}

TEST(Utils, IncrementCounter) {
    /// [Increment counter usage]
    Vec<int, 3> limit = {1, 2, 3};
    Vec<int, 3> cnt; cnt.setZero();

    std::vector<Vec<int, 3>> result;
    do {
        result.push_back(cnt);
    } while (incrementCounter(cnt, limit));
    /// [Increment counter usage]

    std::vector<Vec<int, 3>> expected = {
            {0, 0, 0}, {0, 0, 1}, {0, 0, 2}, {0, 1, 0}, {0, 1, 1}, {0, 1, 2}};
    EXPECT_EQ(expected, result);
}

TEST(Utils, IncrementCounterBothLimits) {
    Vec<int, 3> low = {-3, 2, 0};
    Vec<int, 3> high = {-2, 4, 3};
    Vec<int, 3> cnt = low;

    std::vector<Vec<int, 3>> result;
    do {
        result.push_back(cnt);
    } while (incrementCounter(cnt, low, high));

    std::vector<Vec<int, 3>> expected = {
            {-3, 2, 0}, {-3, 2, 1}, {-3, 2, 2}, {-3, 3, 0}, {-3, 3, 1}, {-3, 3, 2}};
    EXPECT_EQ(expected, result);
}

TEST(Utils, IncrementCounterCyclic) {
    /// [Cyclic counter usage]
    Vec<int, 3> low = {1, 0, 2};
    Vec<int, 3> high = {0, 2, 1};
    Vec<int, 3> size = {3, 3, 3};
    Vec<int, 3> cnt = low;

    std::vector<Vec<int, 3>> result;
    do {
        result.push_back(cnt);
    } while (incrementCyclicCounter(cnt, low, high, size));
    /// [Cyclic counter usage]

    std::vector<Vec<int, 3>> expected = {
            {1, 0, 2}, {1, 0, 0}, {1, 1, 2}, {1, 1, 0}, {2, 0, 2}, {2, 0, 0}, {2, 1, 2}, {2, 1, 0}};
    EXPECT_EQ(expected, result);
}

TEST(Utils, NumLinspace1DBorder) {
    Range<double> expected = {1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2};
    Range<Vec1d> result = linspace(Vec1d(1.0), Vec1d(2.0), 11);
    ASSERT_EQ(expected.size(), result.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_DOUBLE_EQ(expected[i], result[i][0]);;
    }
}

TEST(Utils, NumLinspace1DBorderReversed) {
    Range<double> expected = {2, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 1};
    Range<Vec1d>result = linspace(Vec1d(2.0), {1.0}, 11);
    ASSERT_EQ(expected.size(), result.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_DOUBLE_EQ(expected[i], result[i][0]);;
    }
}

TEST(Utils, NumLinspaceDoubleBorder) {
    Range<double> expected = {1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2};
    Range<double> resultint = linspace(1.0, 2.0, 11);
    ASSERT_EQ(expected.size(), resultint.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_DOUBLE_EQ(expected[i], resultint[i]);;
    }
}

TEST(Utils, NumLinspace1DNoBorder) {
    Range<double> expected = {1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1};
    Range<Vec1d> result = linspace(Vec1d(2.0), {1.0}, 9, false);
    ASSERT_EQ(expected.size(), result.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_DOUBLE_EQ(expected[i], result[i][0]);;
    }
}

TEST(Utils, NumLinspaceDoubleNoBroder) {
    Range<double> expected = {1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1};
    Range<double> resultint = linspace(2.0, 1.0, 9, false);
    ASSERT_EQ(expected.size(), resultint.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_DOUBLE_EQ(expected[i], resultint[i]);;
    }
}

TEST(Utils, NumLinspace2DBorder) {
    Range<Range<double>> expected = {
            {  1, 3}, {  1, 3.1}, {  1, 3.2}, {  1, 3.3}, {  1, 3.4}, {  1, 3.5},
            {1.1, 3}, {1.1, 3.1}, {1.1, 3.2}, {1.1, 3.3}, {1.1, 3.4}, {1.1, 3.5},
            {1.2, 3}, {1.2, 3.1}, {1.2, 3.2}, {1.2, 3.3}, {1.2, 3.4}, {1.2, 3.5}};
    Range<Vec2d> result = linspace(Vec2d{1, 3}, {1.2, 3.5}, {3, 6});
    ASSERT_EQ(expected.size(), result.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_DOUBLE_EQ(expected[i][0], result[i][0]);;
        EXPECT_DOUBLE_EQ(expected[i][1], result[i][1]);;
    }
}

TEST(Utils, NumLinspace2DNoBorder) {
    Range<Range<double>> expected = {
            {1.1, 3.1}, {1.1, 3.2}, {1.1, 3.3}, {1.1, 3.4}};
    Range<Vec2d> result = linspace(Vec2d{1, 3}, {1.2, 3.5}, {1, 4}, false);
    ASSERT_EQ(expected.size(), result.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_DOUBLE_EQ(expected[i][0], result[i][0]);;
        EXPECT_DOUBLE_EQ(expected[i][1], result[i][1]);;
    }
}

TEST(Utils, NumLinspace2DHalfBorder) {
    Range<Range<double>> expected = {
            {  1, 3.1}, {  1, 3.2}, {  1, 3.3}, {  1, 3.4},
            {1.1, 3.1}, {1.1, 3.2}, {1.1, 3.3}, {1.1, 3.4},
            {1.2, 3.1}, {1.2, 3.2}, {1.2, 3.3}, {1.2, 3.4}};
    Range<Vec2d> result = linspace(Vec2d{1, 3}, {1.2, 3.5}, {3, 4}, {true, false});
    ASSERT_EQ(expected.size(), result.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_DOUBLE_EQ(expected[i][0], result[i][0]);;
        EXPECT_DOUBLE_EQ(expected[i][1], result[i][1]);;
    }
}

TEST(Utils, NumLinspace3DBorder) {
    Range<Range<double>> expected = {
            {   0,   1, 2}, {   0,   1, 2.5}, {   0,   1, 3},
            {   0, 1.5, 2}, {   0, 1.5, 2.5}, {   0, 1.5, 3},
            {   0,   2, 2}, {   0,   2, 2.5}, {   0,   2, 3},
            {-0.5,   1, 2}, {-0.5,   1, 2.5}, {-0.5,   1, 3},
            {-0.5, 1.5, 2}, {-0.5, 1.5, 2.5}, {-0.5, 1.5, 3},
            {-0.5,   2, 2}, {-0.5,   2, 2.5}, {-0.5,   2, 3},
            {  -1,   1, 2}, {  -1,   1, 2.5}, {  -1,   1, 3},
            {  -1, 1.5, 2}, {  -1, 1.5, 2.5}, {  -1, 1.5, 3},
            {  -1,   2, 2}, {  -1,   2, 2.5}, {  -1,   2, 3}};
    Range<Vec3d> result = linspace(Vec3d{0, 1, 2}, {-1, 2, 3}, {3, 3, 3});
    ASSERT_EQ(expected.size(), result.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_DOUBLE_EQ(expected[i][0], result[i][0]);;
        EXPECT_DOUBLE_EQ(expected[i][1], result[i][1]);;
        EXPECT_DOUBLE_EQ(expected[i][2], result[i][2]);;
    }
}

TEST(Utils, NumLinspace3DNoBorder) {
    Range<Range<double>> expected = {{-0.5, 1.5, 2.5}};
    Range<Vec3d> result = linspace(Vec3d{0, 1, 2}, {-1, 2, 3}, {1, 1, 1}, false);
    ASSERT_EQ(expected.size(), result.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_DOUBLE_EQ(expected[i][0], result[i][0]);;
        EXPECT_DOUBLE_EQ(expected[i][1], result[i][1]);;
        EXPECT_DOUBLE_EQ(expected[i][2], result[i][2]);;
    }
}

TEST(Utils, NumLinspace3DHalfBorder) {
    Range<Range<double>> expected = {
            {-0.5,   1, 2.5}, {-0.5, 1.5, 2.5}, {-0.5,   2, 2.5}};
    Range<Vec3d> result = linspace(Vec3d{0, 1, 2}, {-1, 2, 3}, {1, 3, 1}, {false, true, false});
    ASSERT_EQ(expected.size(), result.size());
    for (int i = 0; i < expected.size(); ++i) {
        EXPECT_DOUBLE_EQ(expected[i][0], result[i][0]);;
        EXPECT_DOUBLE_EQ(expected[i][1], result[i][1]);;
        EXPECT_DOUBLE_EQ(expected[i][2], result[i][2]);;
    }
}

TEST(Utils, NumBisectionSolve) {
    auto func = [](double x) -> double {
        return (3*std::sin(x*x + 0.4) + 0.1);};
    double target = 2.35;
    double a = bisection(func, 0.0, 1.0, target, 1e-6, 20);
    EXPECT_NEAR(target, func(a), 1e-5);
}

TEST(Utils, NumBisectionHitMaxIterations) {
    auto func = [](double x) -> double {
        return (3*std::sin(x*x + 0.4) + 0.1);};
    double target = 2.35;
    double a = bisection(func, 0.0, 1.0, target, 1e-6, 5);  // must print error message
    EXPECT_NEAR(target, func(a), 1e-2);
}

}  // namespace mm
