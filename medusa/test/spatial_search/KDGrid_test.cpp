#include <medusa/bits/spatial_search/KDGrid_fwd.hpp>

#include "gtest/gtest.h"

#include <medusa/bits/types/Vec.hpp>

namespace mm {

TEST(KDGrid, 2D) {
    KDGrid<Vec2d> search(0, 1, 0.1);

    search.insert(Vec2d(0, 0));
    search.insert(Vec2d(0.05, 0.19));

    EXPECT_TRUE(search.existsPointInSphere(0.0, 1.0));
    EXPECT_FALSE(search.existsPointInSphere(1.0, 1.0));
}

TEST(KDGrid, 3D) {
    KDGrid<Vec2d> search(0, 1, 0.1);

    EXPECT_FALSE(search.existsPointInSphere(0.0, 1.0));

    search.insert(Vec2d(0, 0));
    EXPECT_TRUE(search.existsPointInSphere(0.0, 1.0));
    EXPECT_FALSE(search.existsPointInSphere(1.0, 1.0));
}

TEST(KDGrid, Usage) {
    /// [KDGrid usage example]
    KDGrid<Vec2d> search(0, 1, 0.1);
    search.insert(Vec2d(0, 0));
    if (search.existsPointInSphere(0.0, 1.0)) {
        // do sth ...
    }
    std::cout << search << std::endl;
    /// [KDGrid usage example]
}

}  // namespace mm
