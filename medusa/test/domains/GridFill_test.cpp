#include <medusa/bits/domains/BoxShape_fwd.hpp>
#include <medusa/bits/domains/BallShape_fwd.hpp>
#include <medusa/bits/io/HDF.hpp>
#include "medusa/bits/domains/GridFill_fwd.hpp"
#include "medusa/bits/domains/ShapeDifference.hpp"

#include "gtest/gtest.h"

namespace mm {

TEST(DomainEngines, GridFill2d) {
    /// [GridFill usage example]
    auto sh = BoxShape<Vec2d>(-1, 1) - BallShape<Vec2d>(0, 0.5);
    DomainDiscretization<Vec2d> domain(sh);
    domain.fill(GridFill<Vec2d>(-1, 1), 0.5);
    /// [GridFill usage example]

    auto pos = domain.positions();
    std::sort(pos.begin(), pos.end());
    Range<Vec2d> expected;
    std::vector<double> cs = {-1, -0.5, 0, 0.5, 1};
    for (double x : cs) {
        for (double y : cs) {
            if (x == 0 && y == 0) continue;
            expected.emplace_back(x, y);
        }
    }
    EXPECT_EQ(expected, pos);
}

TEST(DomainEngines, GridFill3d) {
    auto sh = BoxShape<Vec3d>(-1, 1) - BallShape<Vec3d>(0, 0.5);
    DomainDiscretization<Vec3d> domain(sh);
    domain.addInternalNode({0, 0, 1}, 2);  // Test that this node is not doubled.
    domain.fill(GridFill<Vec3d>(-1, 1), 0.5);

    auto pos = domain.positions();
    std::sort(pos.begin(), pos.end());
    Range<Vec3d> expected;
    std::vector<double> cs = {-1, -0.5, 0, 0.5, 1};
    for (double x : cs) {
        for (double y : cs) {
            for (double z : cs) {
                if (x == 0 && y == 0 && z == 0) continue;
                expected.emplace_back(x, y, z);
            }
        }
    }
    EXPECT_EQ(expected, pos);
}


}  // namespace mm
