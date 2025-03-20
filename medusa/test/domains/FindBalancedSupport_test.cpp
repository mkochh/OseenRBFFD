#include <medusa/bits/domains/FindBalancedSupport.hpp>
#include <medusa/bits/types/Vec_fwd.hpp>
#include <medusa/bits/domains/BoxShape.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(DomainEngines, BalancedSupport) {
    Range<Vec2d> points = {
            // find for
            {0,    0},
            {-0.5, 0.5},
            {1,    1},
            // search among
            {1,    1.1},   // 3
            {1,    -1},    // 4
            {-0.5, 0.5},   // 5
            {-2,   2},    // 6
            {-1.1, -1.1},  // 7
            {-3., 3.},  // 8
            // extra & force self
            {0, 0},  // 9
    };
    DomainDiscretization<Vec2d> domain(BoxShape<Vec2d>({-3, -3}, {3, 3}));
    Vec2d shift(0.4, -0.8);  // not centred at 0, 0
    for (const auto& p : points) domain.addInternalNode(p+shift, 1);
    int min_ss = 4, max_ss = 6;
    /// [FindBalancedSupport usage example]
    FindBalancedSupport find_support(min_ss, max_ss);
    find_support.forNodes({0, 1, 2}).searchAmong({3, 4, 5, 6, 7, 8});
    domain.findSupport(find_support);
    /// [FindBalancedSupport usage example]

    for (int i = 0; i < domain.size(); ++i) {
        if (i < 3) {  // test for nodes
            EXPECT_GE(domain.support(i).size(), min_ss);
            EXPECT_LE(domain.support(i).size(), max_ss);
            for (int j : domain.support(i)) {  // test search among
                EXPECT_GE(j, 3);
                EXPECT_LE(j, 8);
            }
        } else {
            EXPECT_TRUE(domain.support(i).empty());
        }
    }
    Range<int> s;
    s = {5, 4, 3, 7};
    EXPECT_EQ(s, domain.support(0));
    s = {5, 3, 7, 4, 6 };
    EXPECT_EQ(s, domain.support(1));
    s = {3, 5, 4, 7, 6, 8};
    EXPECT_EQ(s, domain.support(2));

    find_support.forNodes({9}).searchAmong({3, 4, 5, 6, 7, 8}).forceSelf();
    domain.findSupport(find_support);
    s = {9, 5, 4, 3, 7};
    EXPECT_EQ(s, domain.support(9));
}

TEST(DomainEngines, BalancedSupportBoundary) {
    Range<std::pair<Vec2d, Vec2d>> points_bnd = {  // point and normal
            {{0, 0}, {0, 1}},
            {{0.5, 0}, {0, 1}},
            {{1, 0}, {0, 1}},
            {{-1.1, 0}, {0, 1}},
            {{-1.5, 0}, {0, 1}},
    };
    Range<Vec2d> points = {{-0.5, -0.5}, {0, -0.6}, {0.5, -0.6}};
    DomainDiscretization<Vec2d> domain(BoxShape<Vec2d>({-3, -3}, {3, 3}));
    Vec2d shift(0.4, -0.8);  // not centred at 0, 0
    for (const auto& p : points_bnd) domain.addBoundaryNode(p.first+shift, -1, p.second);
    for (const auto& p : points) domain.addInternalNode(p+shift, 1);
    int min_ss = 4, max_ss = 8;
    FindBalancedSupport find_support(min_ss, max_ss);
    find_support.forNodes({0});
    domain.findSupport(find_support);

    Range<int> s = {0, 1, 6, 5, 7, 2, 3};
    EXPECT_EQ(s, domain.support(0));
}

TEST(DomainEngines, BalancedSupportGhost) {
    DomainDiscretization<Vec3d> domain(UnknownShape<Vec3d>{});
    domain.addBoundaryNode({0, 0, 0}, -1, {-1, 0, 0});
    domain.addGhostNodes(0.1);
    FindBalancedSupport find_support(1, 2);
    domain.findSupport(find_support);
    Range<Range<int>> expected = {{0, 1}, {}};
    EXPECT_EQ(expected, domain.supports());
}

}  // namespace mm
