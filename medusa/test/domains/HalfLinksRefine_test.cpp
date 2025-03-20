#include "medusa/bits/domains/HalfLinksRefine.hpp"

#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/domains/DomainDiscretization.hpp>
#include <medusa/bits/domains/FindClosest.hpp>

#include "gtest/gtest.h"

namespace mm {


TEST(DomainEngines, Refine1d) {
    BoxShape<Vec1d> box(0, 1);
    auto domain = box.discretizeWithStep(0.1);
    int N = domain.size();
    domain.findSupport(FindClosest(3));
    auto region = domain.positions().filter([](const Vec1d &v) { return v[0] < 0.45; });
    HalfLinksRefine refine; refine.region(region);
    auto new_points = refine(domain);
    Range<Vec1d> expected = {0.05, 0.15, 0.25, 0.35, 0.45};
    ASSERT_EQ(expected.size(), new_points.size());
    for (int i = 0; i < expected.size(); ++i) {
        ASSERT_EQ(N + i, new_points[i]);
        EXPECT_DOUBLE_EQ(expected[i][0], domain.pos(new_points[i], 0));
    }
}

template <typename vec_t>
void check_approx_array_equal(
        const Range<vec_t>& expected_pos,
        const Range<int>& expected_types,
        const Range<vec_t>& actual_pos,
        const Range<int>& actual_types,
        double tol) {
//    prn(actual_pos)
//    prn(actual_types)
    ASSERT_EQ(expected_pos.size(), actual_pos.size());
    ASSERT_EQ(expected_types.size(), actual_types.size());
    ASSERT_EQ(expected_pos.size(), expected_types.size());
    int n = expected_pos.size();
    std::vector<bool> seen(n, false);
    for (int i = 0; i < n; ++i) {
        bool found = false;
        for (int j = 0; j < n; ++j) {
            if (!seen[j] && (expected_pos[j] - actual_pos[i]).norm() < tol) {
                seen[j] = true;
                EXPECT_EQ(expected_types[j], actual_types[i]);
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found) << "Element " << actual_pos[i] << " was not found in "
                           << expected_pos << "\n";
    }
    bool all = true;
    for (int i = 0; i < n; ++i) {
        if (!seen[i]) {
            all = false;
            break;
        }
    }
    EXPECT_TRUE(all);
}


TEST(DomainEngines, Refine2dCorner) {
    /// [Half links usage example]
    BoxShape<Vec2d> box(0, 1);
    auto domain = box.discretizeWithStep(0.2);
    domain.findSupport(FindClosest(9));
    int N = domain.size();

    auto region = domain.positions().filter(
            [](const Vec2d &v) { return v[0] < 0.35 && v[1] < 0.35; });
    HalfLinksRefine refine; refine.region(region);
    auto new_points = refine(domain);
    /// [Half links usage example]

    EXPECT_EQ(N + new_points.size(), domain.size());
    Range<Vec2d> expected = {{0.0, 0.1}, {0.0, 0.3}, {0.1, 0.0}, {0.1, 0.1},
                             {0.1, 0.2}, {0.1, 0.3}, {0.2, 0.1}, {0.2, 0.3},
                             {0.3, 0.0}, {0.3, 0.1}, {0.3, 0.2}, {0.3, 0.3}};
    Range<int> expected_types = {-1, -1, -1, 1, 1, 1, 1, 1, -3, 1, 1, 1};
    check_approx_array_equal(expected, expected_types, domain.positions()[new_points].asRange(),
            domain.types()[new_points].asRange(), 1e-10);
}

TEST(DomainEngines, Refine2dLeftSide) {
    BoxShape<Vec2d> box(0, 1);
    auto domain = box.discretizeWithStep(0.5);
    domain.findSupport(FindClosest(9));
    auto region = domain.positions().filter([](const Vec2d &v) { return v[0] < 0.4; });
    HalfLinksRefine refine; refine.region(region);
    auto new_points = refine(domain);
    Range<Vec2d> expected = {{0.0, 0.25}, {0.0, 0.75},  {0.25, 0.0}, {0.25, 0.25},
                             {0.25, 0.5}, {0.25, 0.75}, {0.25, 1.0}};
    Range<int> expected_types = {-1, -1, -1, 1, 1, 1, -4};
    check_approx_array_equal(expected, expected_types, domain.positions()[new_points].asRange(),
            domain.types()[new_points].asRange(), 1e-10);
}

}  // namespace mm
