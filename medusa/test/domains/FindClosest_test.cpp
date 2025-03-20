#include <medusa/bits/types/Vec.hpp>
#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/utils/stdtypesutils.hpp>
#include <map>
#include "medusa/bits/domains/FindClosest.hpp"

#include "gtest/gtest.h"

namespace mm {

TEST(DomainEngines, FindClosestAll) {
    BoxShape<Vec2d> box({0, 0}, {1, 1});
    DomainDiscretization<Vec2d> d = box.discretizeWithStep(0.1);
    int N = d.size();
    Range<int> expected(N);
    for (int i = 0; i < N; ++i) expected[i] = i;
    d.findSupport(FindClosest(N));  // should get everything
    for (auto x : d.supports()) {
        sort(x);
        ASSERT_EQ(expected.size(), x.size());
        EXPECT_EQ(expected, x);
    }
}

// find in radius
/*
// see img/test_case_domain_support.png
TEST(DomainEngines, FindClosestNoLimit) {
    BoxShape<Vec2d> box({0, 0}, {3, 3});
    DomainDiscretization<Vec2d> d = box.discretizeWithStep(1.0);
    d.findSupport(2, 10);  // should not be limiting
    int num = d.size();
    std::map<Vec2d, Range<Vec2d>> expected({
       {{0, 0}, {{0, 0}, {0, 1}, {1, 0}, {1, 1}}},
       {{1, 0}, {{0, 0}, {0, 1}, {1, 0}, {1, 1}, {2, 0}, {2, 1}}},
       {{2, 0}, {{1, 0}, {1, 1}, {2, 0}, {2, 1}, {3, 0}, {3, 1}}},
       {{3, 0}, {{2, 0}, {2, 1}, {3, 0}, {3, 1}}},
       {{0, 1}, {{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}}},
       {{1, 1}, {{0, 0}, {0, 1}, {0, 2}, {1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}}},
       {{2, 1}, {{1, 0}, {1, 1}, {1, 2}, {2, 0}, {2, 1}, {2, 2}, {3, 0}, {3, 1}, {3, 2}}},
       {{3, 1}, {{2, 0}, {2, 1}, {2, 2}, {3, 0}, {3, 1}, {3, 2}}},
       {{0, 2}, {{0, 1}, {0, 2}, {0, 3}, {1, 1}, {1, 2}, {1, 3}}},
       {{1, 2}, {{0, 1}, {0, 2}, {0, 3}, {1, 1}, {1, 2}, {1, 3}, {2, 1}, {2, 2}, {2, 3}}},
       {{2, 2}, {{1, 1}, {1, 2}, {1, 3}, {2, 1}, {2, 2}, {2, 3}, {3, 1}, {3, 2}, {3, 3}}},
       {{3, 2}, {{2, 1}, {2, 2}, {2, 3}, {3, 1}, {3, 2}, {3, 3}}},
       {{0, 3}, {{0, 2}, {0, 3}, {1, 2}, {1, 3}}},
       {{1, 3}, {{0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 2}, {2, 3}}},
       {{2, 3}, {{1, 2}, {1, 3}, {2, 2}, {2, 3}, {3, 2}, {3, 3}}},
       {{3, 3}, {{2, 2}, {2, 3}, {3, 2}, {3, 3}}},
    });
    for (int i = 0; i < num; ++i) {
        Vec2d cur = d.pos(i);
        Range<Vec2d> a = d.supportNodes(i);
        sort(a.begin(), a.end());
        EXPECT_EQ(expected[d.pos(i)], a);
    }
} */


// see img/test_case_domain_support.png
TEST(DomainEngines, FindClosest2d) {
    BoxShape<Vec2d> box({0, 0}, {3, 3});
    DomainDiscretization<Vec2d> d = box.discretizeWithStep(1.0);
    int ss = 5;
    d.findSupport(FindClosest(ss));
    std::map<Vec2d, Range<Vec2d>> expected({
       {{1, 1}, {{0, 1}, {1, 0}, {1, 1}, {1, 2}, {2, 1}}},
       {{2, 1}, {{1, 1}, {2, 0}, {2, 1}, {2, 2}, {3, 1}}},
       {{1, 2}, {{0, 2}, {1, 1}, {1, 2}, {1, 3}, {2, 2}}},
       {{2, 2}, {{1, 2}, {2, 1}, {2, 2}, {2, 3}, {3, 2}}},
    });
    for (const auto& x : d.supports()) {
        EXPECT_LE(x.size(), ss);
    }
    for (int i : d.interior()) {  // internals are well defined
        Vec2d cur = d.pos(i);
        Range<Vec2d> a = d.supportNodes(i);
        sort(a.begin(), a.end());
        EXPECT_EQ(expected[cur], a);
    }
}


class FindClosestTest : public ::testing::Test {
  public:
    DomainDiscretization<Vec1d> domain;
    Range<Vec1d> pos;
    Range<int> types;

    FindClosestTest() : domain(BoxShape<Vec1d>{0, 1}) {}

  protected:
    void SetUp() override {
        pos = {1., 2., 3., 4., 5., 6., 7.};
        types = {-1, 2, -1, 2, 2, -1, -1};

        for (int i = 0; i < pos.size(); ++i) {
            if (types[i] > 0) domain.addInternalNode(pos[i], types[i]);
            else domain.addBoundaryNode(pos[i], types[i], pos[i]);
        }
    }
};

TEST_F(FindClosestTest, Ghost) {
    auto gh = domain.addGhostNodes(0.001);
    domain.findSupport(FindClosest(2));
    Range<Range<int>> expected = {{0, 7}, {1, 7}, {2, 8}, {3, 8}, {4, 3}, {5, 9}, {6, 10},
                                  {}, {}, {}, {}};
    EXPECT_EQ(expected, domain.supports());
}

TEST_F(FindClosestTest, FindClosestForWhich) {
    FindClosest f(1); f.forNodes(domain.boundary());
    domain.findSupport(f);
    Range<Range<int>> expected = {{0}, {}, {2}, {}, {}, {5}, {6}};
    EXPECT_EQ(expected, domain.supports());
}

TEST_F(FindClosestTest, FindClosestSearchAmong) {
    FindClosest f(1); f.searchAmong(domain.interior());
    domain.findSupport(f);
    Range<Range<int>> expected = {{1}, {1}, {1}, {3}, {4}, {4}, {4}};
    EXPECT_EQ(expected, domain.supports());
}

TEST_F(FindClosestTest, FindClosestForceSelf) {
    /// [FindClosest usage example]
    // Find the closest interior node for all nodes.
    FindClosest f(1); f.searchAmong(domain.interior()).forceSelf(false);
    domain.findSupport(f);  // domain is some DomainDiscretization
    /// [FindClosest usage example]
    Range<Range<int>> expected = {{1}, {1}, {1}, {3}, {4}, {4}, {4}};
    EXPECT_EQ(expected, domain.supports());
    f.forceSelf();
    domain.findSupport(f);
    expected = {{0}, {1}, {2}, {3}, {4}, {5}, {6}};
    EXPECT_EQ(expected, domain.supports());
    f.numClosest(2);
    domain.findSupport(f);
    expected = {{0, 1}, {1, 3}, {2, 1}, {3, 4}, {4, 3}, {5, 4}, {6, 4}};
    EXPECT_EQ(expected, domain.supports());
}

TEST(DomainEngines, FindClosestDeath) {
    BoxShape<Vec2d> box({0, 0}, {1, 1});
    DomainDiscretization<Vec2d> d = box.discretizeWithStep(0.1);

    FindClosest f(-2);
    EXPECT_DEATH(d.findSupport(f), "Support size must be greater than 0");
    f.numClosest(122);
    EXPECT_DEATH(d.findSupport(f), "Support size \\(122) cannot exceed number of points that we "
                                   "are searching among \\(121).");
    f.numClosest(5).forNodes({1000});
    EXPECT_DEATH(d.findSupport(f), "Index 1000 out of range \\[0, 121) in forNodes.");
    f.numClosest(1).forNodes({}).searchAmong({5000});
    EXPECT_DEATH(d.findSupport(f), "Index 5000 out of range \\[0, 121) in searchAmong.");
}

}  // namespace mm
