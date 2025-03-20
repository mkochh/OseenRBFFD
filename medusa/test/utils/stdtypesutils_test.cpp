#include <medusa/bits/utils/stdtypesutils.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(Utils, StdTypesSort) {
    std::vector<int> a = {3, 1, 2};
    std::vector<int> b = sorted(a);
    EXPECT_EQ(3, a[0]);
    EXPECT_EQ(1, a[1]);
    EXPECT_EQ(2, a[2]);

    EXPECT_EQ(1, b[0]);
    EXPECT_EQ(2, b[1]);
    EXPECT_EQ(3, b[2]);

    sort(a);
    EXPECT_EQ(b, a);

    b = sorted(a, [](int i, int j) { return j < i; });
    EXPECT_EQ(1, a[0]);
    EXPECT_EQ(2, a[1]);
    EXPECT_EQ(3, a[2]);

    EXPECT_EQ(3, b[0]);
    EXPECT_EQ(2, b[1]);
    EXPECT_EQ(1, b[2]);

    sort(a, [](int i, int j) { return j < i; });
    EXPECT_EQ(b, a);
}

TEST(Utils, StdTypesPad) {
    std::vector<std::vector<int>> a = {{2, 1}, {4}, {}, {1, 2, 3}};
    std::vector<std::vector<int>> expected = {{2, 1, -1}, {4, -1, -1}, {-1, -1, -1}, {1, 2, 3}};
    EXPECT_EQ(expected, pad(a, -1));
    a = {{1, 2}, {3, 4}, {5, 6}};
    EXPECT_EQ(a, pad(a, -1));
}

TEST(Utils, StdTypesTextSplit) {
    EXPECT_EQ(std::vector<std::string>({"a", "b", "c"}), split("a,b,c", ','));
    EXPECT_EQ(std::vector<std::string>({"a", "b", "c"}), split("a, b, c", ", "));
    EXPECT_EQ(std::vector<std::string>({"", "", "", " ", ""}), split("aaa a", "a"));
    EXPECT_EQ(std::vector<std::string>({"a", "b", ""}), split("a,b,", ","));
    EXPECT_EQ(std::vector<std::string>({"a", "b", ""}), split("a,,b,,", ",,"));
}

TEST(Utils, StdTypesTextJoin) {
    EXPECT_EQ("ababa", join({"a", "a", "a"}, 'b'));
    EXPECT_EQ("ababa", join({"a", "a", "a"}, "b"));
    EXPECT_EQ("abxabxa", join({"a", "a", "a"}, "bx"));
    EXPECT_EQ("aa", join({"a", "", "a"}, ""));
    EXPECT_EQ("a123123123123a123", join({"a", "", "", "", "a", ""}, "123"));
}

TEST(Utils, PrintTuple) {
    std::stringstream ss;
    std::tuple<int, double, std::string> t(2, 4.5, "asdf");
    ss << t;
    EXPECT_EQ("(2, 4.5, asdf)", ss.str());

    ss.str({});
    ss.clear();
    std::tuple<> empty;
    ss << empty;
    EXPECT_EQ("()", ss.str());
}

}  // namespace mm
