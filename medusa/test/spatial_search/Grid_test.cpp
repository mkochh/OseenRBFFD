#include <medusa/bits/spatial_search/Grid.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(Grid, 3D) {
    Grid<int, 3> grid({2, 3, 4}, -1);

    int c = 0;
    for (int i = 0; i < grid.size(0); ++i) {
        for (int j = 0; j < grid.size(1); ++j) {
            for (int k = 0; k < grid.size(2); ++k) {
                grid({i, j, k}) = c++;
            }
        }
    }

    // test setting, getting and row major storage
    for (int i = 0; i < grid.size(0); ++i) {
        for (int j = 0; j < grid.size(1); ++j) {
            for (int k = 0; k < grid.size(2); ++k) {
                ASSERT_EQ(grid.linearIndex({i, j, k}), grid({i, j, k}));
            }
        }
    }
}

TEST(Grid, 1D) {
    Grid<double, 1> grid({10});
    grid({0}) = -3;
    grid({3}) = 2;
    EXPECT_EQ(-3, grid({0}));
    EXPECT_EQ(0, grid({1}));
    EXPECT_EQ(0, grid({2}));
    EXPECT_EQ(2, grid({3}));
    EXPECT_EQ(0, grid({4}));
}

TEST(Grid, Index3D) {
    Grid<double, 3> grid({3, 4, 5});

    EXPECT_EQ(0, grid.linearIndex({0, 0, 0}));
    EXPECT_EQ(4, grid.linearIndex({0, 0, 4}));
    EXPECT_EQ(19, grid.linearIndex({0, 3, 4}));

    EXPECT_EQ((std::array<int, 3>({{0, 0, 0}})), grid.multiIndex(0));
    EXPECT_EQ((std::array<int, 3>({{0, 0, 4}})), grid.multiIndex(4));
    EXPECT_EQ((std::array<int, 3>({{0, 3, 4}})), grid.multiIndex(19));
}

TEST(Grid, Usage_4D) {
    /// [Grid usage example]
    Grid<double, 4> grid({2, 2, 2, 2});
    grid({0, 1, 0, 0}) = -3;
    grid({0, 0, 1, 0}) = 2;
    EXPECT_EQ(2*2*2*2, grid.size());
    EXPECT_EQ((std::array<int, 4>({2, 2, 2, 2})), grid.sizes());
    EXPECT_EQ(2, grid.size(2));
    EXPECT_EQ(-3, grid({0, 1, 0, 0}));
    EXPECT_EQ(2, grid({0, 0, 1, 0}));
    std::cout << grid << std::endl;
    /// [Grid usage example]
}

}  // namespace mm
