#include <medusa/bits/types/Range.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(Types, RangeConstructAssign) {
    Range<double> a({1, 2, 3.6});
    ASSERT_EQ(3, a.size());
    EXPECT_EQ(1, a[0]);
    EXPECT_EQ(2, a[1]);
    EXPECT_EQ(3.6, a[2]);
    a = {2.3, -4.1};
    ASSERT_EQ(2, a.size());
    EXPECT_EQ(2.3, a[0]);
    EXPECT_EQ(-4.1, a[1]);
    a = {3, 4, 5, 6, 7.2};
    Range<double> b(1, 3);
    b = a[{1, 2, 4}].asRange();
    ASSERT_EQ(3, b.size());
    EXPECT_EQ(4, b[0]);
    EXPECT_EQ(5, b[1]);
    EXPECT_EQ(7.2, b[2]);
}

TEST(Types, RangeCompare) {
    Range<double> a = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    Range<int> result = a > 3.2;
    EXPECT_EQ(Range<int>({3, 4, 5, 6, 7}), result);
    result = a < 3.2;
    EXPECT_EQ(Range<int>({0, 1, 2}), result);
    result = a == 3.2;
    EXPECT_EQ(Range<int>(), result);
    EXPECT_EQ(a, a);

    Range<double> b(std::vector<double>({1, 2, 3}));
    ASSERT_EQ(3, b.size());
    EXPECT_EQ(1, b[0]);
    EXPECT_EQ(2, b[1]);
    EXPECT_EQ(3, b[2]);

    b = std::vector<double>({1.2, 2.3});
    ASSERT_EQ(2, b.size());
    EXPECT_DOUBLE_EQ(1.2, b[0]);
    EXPECT_DOUBLE_EQ(2.3, b[1]);

    Range<int> c = {1, 2, 3, 1, 2, 3};
    Range<int> d = {1, 2, 3};
    EXPECT_EQ((c[{0, 1, 2}]).asRange(), d);
    EXPECT_NE((c[{1, 2, 0}]).asRange(), d);
    EXPECT_EQ((c[{0, 1, 2}]).asRange(), (c[{3, 4, 5}]).asRange());
    EXPECT_EQ(d, (c[{0, 1, 2}]).asRange());
}

TEST(Types, RangeAccess) {
    Range<double> d_range(5, 1);
    d_range = static_cast<Range<double>>(d_range[{0, 2}]);
    ASSERT_EQ(2, d_range.size());
    EXPECT_EQ(1, d_range[0]);
    EXPECT_EQ(1, d_range[1]);
}

TEST(Types, RangeAppendJoin) {
    Range<double> a = {1.0, -2.2, 3.0}, b = {2.0, 3.0, 4.0};
    Range<double> c(b);
    Range<double> d(a);
    a.append(b);
    EXPECT_EQ(Range<double>({1.0, -2.2, 3.0, 2.0, 3.0, 4.0}), a);
    EXPECT_EQ(c, b);
    EXPECT_EQ(Range<double>({1.0, -2.2, 3.0, 2.0, 3.0, 4.0}), d.join(c));
}

TEST(Types, RangeRemove) {
    Range<double> a = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    a.remove({4, 3, 2, 3});
    EXPECT_EQ(Range<double>({0, 1, 5, 6, 7, 8}), a);
    a.remove({});
    EXPECT_EQ(Range<double>({0, 1, 5, 6, 7, 8}), a);
    a.remove({0, 1, 2, 3, 4, 5});
    EXPECT_EQ(Range<double>(), a);
}

TEST(Types, RangeFilter) {
    Range<double> aa = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    Range<int> ind = aa.filter([](double v) { return v > 2.3 && v < 6.4; });
    EXPECT_EQ(Range<int>({2, 3, 4, 5}), ind);
}

TEST(Types, RangeMap) {
    Range<std::pair<double, double>> r = {{1, 0}, {1.5, 2}, {0, -0.1}};
    Range<double> n = r.map([] (const std::pair<double, double>& v) {
        return v.first*v.first + v.second*v.second; });
    ASSERT_EQ(3, n.size());
    EXPECT_DOUBLE_EQ(1.0*1.0, n[0]);
    EXPECT_DOUBLE_EQ(2.5*2.5, n[1]);
    EXPECT_DOUBLE_EQ(0.1*0.1, n[2]);
}

TEST(Types, RangeBool) {
    Range<bool> a(30, false);
    a[7] = true;
    EXPECT_TRUE(a[7]);
    EXPECT_FALSE(a[6]);
    EXPECT_FALSE(a[8]);
}

TEST(Types, RangeDeathTest) {
     Range<int> a = {1, 2, 3};
     Range<int> b = {3, 1};
     EXPECT_DEATH((a[{4, 1}]),
                  "Index 4 out of range \\[0, 3) when using multiindexed read-write access.");
     EXPECT_DEATH((a[{4, 1}]) = 5,
                  "Index 4 out of range \\[0, 3) when using multiindexed read-write access.");
     EXPECT_DEATH((static_cast<const Range<int>&>(a)[{4, 1}]),
                  "Index 4 out of range \\[0, 3) when using multiindexed read access.");
     EXPECT_DEATH(static_cast<const Range<int>&>(a)[4],
                  "Index 4 out of range \\[0, 3) when accessing Range for read/write.");
     EXPECT_DEATH(a[4], "Index 4 out of range \\[0, 3) when accessing Range for write.");
}

TEST(Types, RangeViewPrint) {
    std::stringstream ss;
    Range<int> a = {1, 2, 3};
    ss << a[{1, 0, 1}];
    ss << static_cast<const Range<int>&>(a)[{1, 2, 1}];

    EXPECT_EQ("[2, 1, 2][2, 3, 2]", ss.str());
}

TEST(Types, Seq) {
    EXPECT_EQ(Range<int>({0, 1, 2, 3, 4}), Range<int>::seq(5));
    EXPECT_EQ(Range<int>({2, 3, 4}), Range<int>::seq(2, 5));
    EXPECT_EQ(Range<int>({4, 3, 2, 1, 0}), Range<int>::seq(4, -1, -1));
    EXPECT_EQ(Range<int>({4, 2, 0}), Range<int>::seq(4, -1, -2));
    EXPECT_EQ(Range<int>({4, 2}), Range<int>::seq(4, 0, -2));
    EXPECT_EQ(Range<int>({2, 4}), Range<int>::seq(2, 5, 2));
    EXPECT_EQ(Range<int>({2, 4}), Range<int>::seq(2, 6, 2));
    EXPECT_EQ(Range<int>({2, 4, 6}), Range<int>::seq(2, 7, 2));
}

TEST(Types, DISABLED_RangeExamples) {
    /// [Range usage example]
    Range<int> v = Range<int>::seq(100);  // numbers from 0 to 99, inclusive
    v[v.filter([](int x) { return x % 2 == 0; })] = 0;  // set even to 0
    v.remove(v > 80);  // remove greater than 80
    std::cout << v << std::endl;

    Range<int> u(100, 4);
    v.append(u);  // append one range to another
    std::cout << v + u << std::endl;  // join two ranges and return result
    /// [Range usage example]
}

}  // namespace mm
