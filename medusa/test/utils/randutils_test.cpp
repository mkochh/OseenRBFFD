#include <medusa/bits/utils/randutils.hpp>

#include "gtest/gtest.h"
#include <medusa/bits/types/Range.hpp>

namespace mm {

TEST(Utils, RandGetSeed) {
    assert(get_seed() != get_seed());
}

TEST(Utils, RandChoice) {
    // random choice -- deterministic tests
    Range<int> elements = {1, 2, 3};
    Range<double> weights = {0, 3, 0};
    EXPECT_EQ(random_choice(elements, weights), 2);
    weights = {0, 1, 0};
    EXPECT_EQ(random_choice(elements, weights, true), 2);

    elements = {0, 1, 2};
    std::mt19937 generator(42);
    std::vector<int> counts = {0, 0, 0};
    for (int i = 0; i < 100; ++i) {
        int c = random_choice(elements, {1, 1, 1}, false, generator);
        EXPECT_LE(0, c);
        EXPECT_LE(c, 2);
        counts[c]++;
    }
    EXPECT_LE(30, counts[0]);
    EXPECT_LE(counts[0], 37);
    EXPECT_LE(30, counts[1]);
    EXPECT_LE(counts[1], 37);
    EXPECT_LE(30, counts[2]);
    EXPECT_LE(counts[2], 37);
}


}  // namespace mm
