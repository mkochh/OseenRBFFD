#include <medusa/bits/domains/DomainDiscretization.hpp>
#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/domains/BallShape.hpp>
#include <medusa/bits/io/HDF.hpp>
#include <medusa/bits/domains/FindClosest.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(Domains, AddNodes) {
    UnknownShape<Vec2d> sh;
    DomainDiscretization<Vec2d> d(sh);
    d.assert_is_valid();
    int idx = d.addInternalNode({-2, 3}, 3);
    d.assert_is_valid();
    EXPECT_EQ(0, idx);
    idx = d.addBoundaryNode({0, 1}, -2, {-1, 0});
    d.assert_is_valid();
    EXPECT_EQ(1, idx);
}

TEST(Domains, ChangeNodes) {
    DomainDiscretization<Vec1d> domain(UnknownShape<Vec1d>{});
    domain.addBoundaryNode(0, -1, -1);
    domain.addInternalNode(1, 2);
    domain.addInternalNode(2, 1);
    domain.addInternalNode(3, 3);
    domain.addBoundaryNode(4, -5, 1);
//    EXPECT_DEATH(domain.changeToBoundary(-1, 0, 0, 0), "");
//    EXPECT_DEATH(domain.changeToInterior(-1, 0, 0), "");
//    EXPECT_DEATH(domain.changeToBoundary(0, 0, 0, 0), "");
//    EXPECT_DEATH(domain.changeToInterior(1, 0, 0), "");
//    EXPECT_DEATH(domain.changeToBoundary(1, 0, 1, 0), "");
//    EXPECT_DEATH(domain.changeToInterior(0, 0, -1), "");
    std::vector<Vec1d> normals = {{-1}, {1}};
    std::vector<Vec1d> positions = {{0}, {1}};
    domain.changeToBoundary(1, 0, -1, -1);
    domain.assert_is_valid();
    EXPECT_EQ(normals[0], domain.normal(1));
    EXPECT_EQ(positions[0], domain.pos(1));
    EXPECT_EQ(-1, domain.type(1));
    domain.changeToInterior(0, 1, 1);
    domain.assert_is_valid();
    EXPECT_EQ(positions[1], domain.pos(0));
    EXPECT_EQ(1, domain.type(0));
    EXPECT_EQ(normals[0], domain.normal(1));
    EXPECT_EQ(normals[1], domain.normal(4));
}

TEST(Domains, Accessors) {
    /// [UnknownShape usage example]
    UnknownShape<Vec2d> sh;
    DomainDiscretization<Vec2d> d(sh);
    std::vector<Vec2d> pos = {{-2.3, 3.4}, {-2.1, 3.4}, {0, 1.4}, {-1.4, 0.001}};
    std::vector<int> types = {2, -3, 1, -1};
    std::vector<Vec2d> normals = {{-1.0, 0}, {0.0, 1.0}};
    d.addInternalNode(pos[0], types[0]);
    d.addBoundaryNode(pos[1], types[1], normals[0]);
    d.addInternalNode(pos[2], types[2]);
    d.addBoundaryNode(pos[3], types[3], normals[1]);
    /// [UnknownShape usage example]

    EXPECT_EQ(4, d.size());
    EXPECT_EQ(pos, d.positions());
    EXPECT_EQ(pos[2], d.pos(2));
    EXPECT_EQ(pos[2][1], d.pos(2, 1));

    EXPECT_EQ(types, d.types());
    EXPECT_EQ(types[2], d.type(2));

    std::vector<int> expected = {0, 2};
    EXPECT_EQ(expected, d.interior());
    expected = {1, 3};
    EXPECT_EQ(expected, d.boundary());

    EXPECT_EQ(normals[0], d.normal(1));
    EXPECT_EQ(normals[1], d.normal(3));

    expected = {0, 0, 0, 0};
    EXPECT_EQ(expected, d.supportSizes());

    std::vector<std::vector<int>> supports = {{0, 1}, {1, 0}, {2, 0, 1}, {3, 2}};
    for (int i = 0; i < 4; ++i) d.support(i) = supports[i];
    expected = {2, 2, 3, 2};
    EXPECT_EQ(expected, d.supportSizes());
    EXPECT_EQ(2, d.supportSize(0));
    EXPECT_EQ(3, d.supportSize(2));
    EXPECT_EQ(supports[3], d.support(3));
    EXPECT_EQ(supports[2][1], d.support(2, 1));
    EXPECT_EQ(pos[supports[2][1]], d.supportNode(2, 1));
    EXPECT_DOUBLE_EQ((pos[0] - pos[1]).norm(), d.dr(1));
    EXPECT_DOUBLE_EQ((pos[0] - pos[1]).norm(), d.dr(0));
}

TEST(Domains, Reorder) {
    DomainDiscretization<Vec1d> domain(UnknownShape<Vec1d>{});
    domain.addBoundaryNode(-1, -1, -1);
    domain.addInternalNode(2, 2);
    domain.addInternalNode(1, 1);
    domain.addInternalNode(3, 3);
    domain.addBoundaryNode(5, -5, 1);
    domain.findSupport(FindClosest(5));
    auto domain_old = domain;

    indexes_t p = domain.reorderNodes();
    domain.assert_is_valid();

    for (int i = 0; i < 5; ++i) {
        if (i > 0) {
            EXPECT_LT(domain.pos(i-1, 0), domain.pos(i, 0));
        }
        EXPECT_EQ(domain_old.pos(p[i]), domain.pos(i));
        EXPECT_EQ(domain_old.type(p[i]), domain.type(i));
        if (domain.type(i) < 0) {
            EXPECT_EQ(domain_old.normal(p[i]), domain.normal(i));
        }
        EXPECT_EQ(domain_old.supportNodes(p[i]), domain.supportNodes(i));
    }
}

TEST(Domains, ReorderFn) {
    DomainDiscretization<Vec1d> domain(UnknownShape<Vec1d>{});
    domain.addBoundaryNode(-1, -1, -1);
    domain.addInternalNode(2, 2);
    domain.addInternalNode(1, 1);
    domain.addInternalNode(3, 3);
    domain.addBoundaryNode(5, -5, 1);
    domain.findSupport(FindClosest(5));
    auto domain_old = domain;

    indexes_t p = domain.reorderNodes(
            [](const std::pair<Vec1d, int>& a, const std::pair<Vec1d, int>& b) {
                return a.first[0] > b.first[0];
            });
    domain.assert_is_valid();

    for (int i = 0; i < 5; ++i) {
        if (i > 0) {
            EXPECT_GT(domain.pos(i-1, 0), domain.pos(i, 0));
        }
        EXPECT_EQ(domain_old.pos(p[i]), domain.pos(i));
        EXPECT_EQ(domain_old.type(p[i]), domain.type(i));
        if (domain.type(i) < 0) {
            EXPECT_EQ(domain_old.normal(p[i]), domain.normal(i));
        }
        EXPECT_EQ(domain_old.supportNodes(p[i]), domain.supportNodes(i));
    }
}

TEST(Domains, Shuffle) {
    DomainDiscretization<Vec1d> domain(UnknownShape<Vec1d>{});
    domain.addBoundaryNode(-1, -1, -1);
    domain.addInternalNode(2, 2);
    domain.addInternalNode(1, 1);
    domain.addInternalNode(3, 3);
    domain.addBoundaryNode(5, -5, 1);
    domain.findSupport(FindClosest(5));
    auto domain_old = domain;

    indexes_t p = {4, 3, 2, 0, 1};
    domain.shuffleNodes(p);

    domain.assert_is_valid();
    for (int i = 0; i < 5; ++i) {
        EXPECT_EQ(domain_old.pos(p[i]), domain.pos(i));
        EXPECT_EQ(domain_old.type(p[i]), domain.type(i));
        if (domain.type(i) < 0) {
            EXPECT_EQ(domain_old.normal(p[i]), domain.normal(i));
        }
        EXPECT_EQ(domain_old.supportNodes(p[i]), domain.supportNodes(i));
    }
}

TEST(Domains, GhostNodes) {
    DomainDiscretization<Vec1d> domain(UnknownShape<Vec1d>{});
    domain.addBoundaryNode(0, -1, -1);
    domain.addBoundaryNode(1, -1, 1);
    domain.addInternalNode(0.5, 1);
    int n = domain.size();
    Range<int> gh = domain.addGhostNodes(0.5);
    int n2 = domain.size();
    EXPECT_EQ(n+2, n2);
    EXPECT_EQ(indexes_t({3, 4, -1}), gh);
    EXPECT_EQ(-0.5, domain.pos(gh[0], 0));
    EXPECT_EQ(1.5, domain.pos(gh[1], 0));
    auto gh2 = domain.addGhostNodes([](const Vec1d& p) { return (p[0] > 0) ? 1 : -1; }, -2, {1});
    int n3 = domain.size();
    EXPECT_EQ(n3-1, gh2[1]);
    EXPECT_EQ(2, domain.pos(n3-1, 0));
    EXPECT_EQ(-2, domain.type(n3-1));
    EXPECT_EQ(1, domain.normal(n3-1)[0]);
}

TEST(Domains, DiscretizeAdd) {
    DomainDiscretization<Vec2d> rect(BoxShape<Vec2d>({0, 0}, {2, 1}));
    rect.addInternalNode({1, 0.5}, 1);
    rect.addBoundaryNode({1.5, 0}, -1, {1, 0.8});
    rect.addInternalNode({1.5, 0.5}, 1);
    rect.addBoundaryNode({2, 0.5}, -1, {0.7, 1});
    rect.addBoundaryNode({0, 0.5}, -1, {-1, 0});
    rect.assert_is_valid();
    DomainDiscretization<Vec2d> union_test(BallShape<Vec2d>({2, 1}, 1));
    union_test.assert_is_valid();
    union_test.addBoundaryNode({2 - std::sqrt(2) / 2, 1 - std::sqrt(2) / 2}, -1, {1, -1});
    union_test.addInternalNode({1.5, 0.5}, 1);
    union_test.addBoundaryNode({3, 1}, -1, {0.4, 6.5});
    union_test.addInternalNode({2, 1.5}, 1);
    union_test.addBoundaryNode({2, 0}, -1, {0, -1});
    union_test.assert_is_valid();
    union_test += rect;
    union_test.assert_is_valid();
    ShapeUnion<Vec2d> sh = dynamic_cast<const ShapeUnion<Vec2d>&>(union_test.shape());
    (void) sh;  // Try casting sh to assert correct type.
    Range<Vec2d> expected_boundary = {{3, 1}, {1.5, 0}, {0, 0.5}};
    Range<Vec2d> boundary = union_test.positions()[union_test.boundary()];
    EXPECT_EQ(expected_boundary, boundary);
    Range<Vec2d> expected_interior = {{2, 1.5}, {1, 0.5}, {1.5, 0.5}};
    Range<Vec2d> interior = union_test.positions()[union_test.interior()];
    EXPECT_EQ(expected_interior, interior);
}

TEST(Domains, DiscretizeAddDisjoint) {
    DomainDiscretization<Vec2d> d1 = BoxShape<Vec2d>({0, 0}, {1, 1}).discretizeWithStep(0.1);
    DomainDiscretization<Vec2d> d2 = BoxShape<Vec2d>({2, 2}, {3, 3}).discretizeWithStep(0.1);

    int d1int = d1.interior().size();
    int d1bnd = d1.boundary().size();
    d1 += d2;

    EXPECT_EQ(d1int+d2.interior().size(), d1.interior().size());
    EXPECT_EQ(d1bnd+d2.boundary().size(), d1.boundary().size());
}

TEST(Domains, DiscretizeAdd1D) {
    double dx = 0.1;
    BallShape<Vec1d> b1(0.5, 0.5);
    BallShape<Vec1d> b2(0, 0.25);
    auto domain1 = b1.discretizeBoundaryWithStep(dx);
    auto domain2 = b2.discretizeBoundaryWithStep(dx);

    domain1 += domain2;

    Range<Vec1d> expected_pos = {{1}, {-0.25}};
    Range<int> expected_types = {-1, -1};
    EXPECT_EQ(expected_types, domain1.types());
    EXPECT_EQ(expected_pos, domain1.positions());
}

TEST(Domains, DiscretizeAddLarge) {
    double h = 2*PI/60;
    DomainDiscretization<Vec2d> rect = BoxShape<Vec2d>({-1, -1}, {1, 1}).discretizeWithStep(h);
    DomainDiscretization<Vec2d> circ = BallShape<Vec2d>({0, 0}, 0.7).discretizeWithStep(h);
    int rsize = rect.size();
    int csize = circ.size();
    auto cbnd = circ.boundary();
    int cbsize = cbnd.size();
    circ.types() = 5;
    circ.types()[cbnd] = -7;
    rect += circ;
    EXPECT_LT(rect.size(), rsize+csize);
    EXPECT_GT(rect.size(), csize);
    KDTree<Vec2d> tree(rect.positions());
    double dx = 1.0;
    for (int i = 0; i < rect.size(); ++i) {
        dx = std::min(dx, std::sqrt(tree.query(rect.pos(i), 2).second[1]));
    }
    EXPECT_GE(dx, 0.099);  // due to amortisation in circle and box boundary discretizations
    EXPECT_EQ(csize - cbsize, (rect.types() == 5).size());
}

TEST(Domains, DiscretizeSubtract) {
    DomainDiscretization<Vec2d> rect(BoxShape<Vec2d>({0, 0}, {2, 1}));
    rect.addInternalNode({1, 0.5}, 1);
    rect.addBoundaryNode({1.5, 0}, -1, {0, -1});
    rect.addBoundaryNode({1.5, 1}, -1, {0, 0.8});
    rect.addInternalNode({1.5, 0.5}, 1);
    rect.addBoundaryNode({2, 0.5}, -1, {0.7, 1});
    rect.addBoundaryNode({0, 0.5}, -1, {-1, 0});
    rect.assert_is_valid();
    DomainDiscretization<Vec2d> union_test(BallShape<Vec2d>({2, 1}, 1));
    union_test.assert_is_valid();
    union_test.addBoundaryNode({2 - std::sqrt(2) / 2, 1 - std::sqrt(2) / 2}, -1, {1, -1});
    union_test.addInternalNode({1.5, 0.5}, 1);
    union_test.addBoundaryNode({3, 1}, -1, {0.4, 6.5});
    union_test.addInternalNode({2, 1.5}, 1);
    union_test.addBoundaryNode({2, 0}, -1, {0, -1});
    union_test.addBoundaryNode({2, 2}, -1, {0, 1});
    union_test.assert_is_valid();
    union_test -= rect;
    ShapeDifference<Vec2d> sh = dynamic_cast<const ShapeDifference<Vec2d>&>(union_test.shape());
    (void) sh;  // Try casting sh to assert correct type.
    union_test.assert_is_valid();
    Range<Vec2d> expected_boundary = {{3, 1}, {2, 2}, {1.5, 1}, {2, 0.5}};
    Range<Vec2d> boundary = union_test.positions()[union_test.boundary()];
    EXPECT_EQ(expected_boundary, boundary);
    Range<Vec2d> expected_interior = {{2, 1.5}};
    Range<Vec2d> interior = union_test.positions()[union_test.interior()];
    EXPECT_EQ(expected_interior, interior);
}

TEST(Domains, DiscretizeSubtractDisjoint) {
    DomainDiscretization<Vec2d> d1 = BoxShape<Vec2d>({0, 0}, {1, 1}).discretizeWithStep(0.1);
    DomainDiscretization<Vec2d> d2 = BoxShape<Vec2d>({2, 2}, {3, 3}).discretizeWithStep(0.1);

    int d1int = d1.interior().size();
    int d1bnd = d1.boundary().size();
    d1 -= d2;

    EXPECT_EQ(d1int, d1.interior().size());
    EXPECT_EQ(d1bnd, d1.boundary().size());
}

TEST(Domains, DiscretizeSubtractLarge) {
    double h = 2*PI/60;
    DomainDiscretization<Vec2d> rect = BoxShape<Vec2d>({-1, -1}, {1, 1}).discretizeWithStep(h);
    DomainDiscretization<Vec2d> circ = BallShape<Vec2d>({0, 0}, 0.7).discretizeWithStep(h);
    int rsize = rect.size();
    auto cbnd = circ.boundary();
    int cbsize = cbnd.size();
    circ.types() = 5;
    circ.types()[cbnd] = -7;
    rect -= circ;
    EXPECT_LT(rect.size(), rsize);
    KDTree<Vec2d> tree(circ.positions());
    double dx = 1.0;
    for (int i = 0; i < circ.size(); ++i) {
        dx = std::min(dx, std::sqrt(tree.query(circ.pos(i), 2).second[1]));
    }
    EXPECT_GE(dx, 0.099);  // due to amortisation in circle and box boundary discretizations
    EXPECT_EQ(cbsize, (rect.types() == -7).size());
}

TEST(Domains, ProjectPointToBoundary1d) {
    BallShape<Vec1d> circ({2.2}, 0.7);
    auto discr = circ.discretizeWithStep(0.1);
    Range<Vec1d> positions = discr.positions();

    bool success;
    Vec1d added_point;
    Vec1d hint(1.5);
    std::tie(success, added_point) = discr.shape().projectPointToBoundary(hint, {-1});
    EXPECT_FALSE(success);
    EXPECT_EQ(positions, discr.positions());
}

TEST(Domains, ProjectPointToBoundary2d) {
    const double precision = 1e-5;
    /// [Add to boundary]
    Vec2d center(2.2, 1.3);
    BallShape<Vec2d> circ(center, 0.7);
    auto discr = circ.discretizeBoundaryWithStep(0.1);

    bool success;
    Vec2d added_point;
    Vec2d hint = {1.5, 0.8};
    Vec2d normal = (hint-center).normalized();
    std::tie(success, added_point) = discr.shape().projectPointToBoundary(hint, normal);
    /// [Add to boundary]

    EXPECT_TRUE(success);
    EXPECT_NEAR(1.6303909089, added_point[0], precision);
    EXPECT_NEAR(0.8931363635, added_point[1], precision);

    hint = {2.8, 1.6};
    normal = (hint-center).normalized();
    std::tie(success, added_point) = discr.shape().projectPointToBoundary(hint, normal);
    EXPECT_TRUE(success);
    EXPECT_NEAR(2.826098822, added_point[0], precision);
    EXPECT_NEAR(1.613049411, added_point[1], precision);

    hint = {-2.8, 1.6};
    normal = {0, 1};
    std::tie(success, added_point) = discr.shape().projectPointToBoundary(hint, normal);
    EXPECT_FALSE(success);
}

TEST(Domains, ProjectPointToBoundary3d) {
    BallShape<Vec3d> circ({0, 0, 0}, 1);
    auto discr = circ.discretizeBoundaryWithStep(0.1);
    double precision = 1e-2;

    bool success;
    Vec3d added_point;

    Vec3d hint = {0.8, 0.01, -0.01};
    std::tie(success, added_point) = discr.shape().projectPointToBoundary(
            hint, hint.normalized());
    EXPECT_TRUE(success);
    EXPECT_NEAR(0.99947, added_point[0], precision);
    EXPECT_NEAR(0.0117135, added_point[1], precision);
    EXPECT_NEAR(-0.0105922, added_point[2], precision);
    EXPECT_NEAR(added_point.squaredNorm(), 1, precision);
}

TEST(Domains, LoadTest) {
    /// [load usage example]
    // Load domain from file `test_domain_load.h5`. It is stored in the `domain` group.
    // Shape of `d` is UnknownShape.
    HDF file("test/testdata/test_domain_load.h5", HDF::READONLY);
    auto d = DomainDiscretization<Vec2d>::load(file, "domain");
    /// [load usage example]

    BoxShape<Vec2d> b(0.0, 1.0);
    DomainDiscretization<Vec2d> expected = b.discretizeWithStep(0.1);
    EXPECT_EQ(expected.positions(), d.positions());
    EXPECT_EQ(expected.types(), d.types());
    EXPECT_EQ(expected.bmap(), d.bmap());
    EXPECT_EQ(expected.normals(), d.normals());
}

TEST(Domains, DiscreteContainsTest) {
    /// [discreteContains usage example]
    BoxShape<Vec3d> box({0, 0, 0}, {2, 2, 2});
    BallShape<Vec3d> circ({0, 0, 0}, 1);
    auto shape = box - circ;
    auto discr = shape.discretizeBoundaryWithStep(0.1);

    KDTree<Vec3d> search;
    discr.makeDiscreteContainsStructure(search);
    EXPECT_TRUE(discr.discreteContains({1.2, 1.2, 1.2}, search));
    EXPECT_FALSE(discr.discreteContains({0, 0, 0}, search));
    // discreteContains might fail near the shape boundary
    EXPECT_FALSE(discr.discreteContains({2, 2, 2}, search));
    /// [discreteContains usage example]

    Range<Vec3d> points({{1.5, 1.5, 1.5}, {1.9, 1.9, 1.9}, {3, 3, 3}, {2.1, 2.1, 2.1},
                         {0.1, 1.7, 0.6}, {1.8, 0.2, 0.1}, {3.1, 16.0, 0.0}});
    for (const auto& p : points) {
        EXPECT_EQ(discr.discreteContains(p, search), discr.contains(p, search));
    }
}


DomainDiscretization<Vec2d> test_bool_op(bool with_interior, bool add) {
    double h = 0.25;
    BallShape<Vec2d> sA({0, 0}, 2);
    BallShape<Vec2d> sB({3, 0}, 2);
    auto A = with_interior ? sA.discretizeWithStep(h) : sA.discretizeBoundaryWithStep(h);
    auto B = with_interior ? sB.discretizeWithStep(h) : sB.discretizeBoundaryWithStep(h);
    // See plot_domain_bool_op.m for plotting.
//    prn("apos", A.positions())
//    prn("at", A.types())
//    prn("bpos", B.positions())
//    prn("bt", B.types())
    if (add) {
        A += B;
    } else {
        A -= B;
    }
//    prn("abpos", A.positions())
//    prn("abt", A.types())
    return A;
}

TEST(Domains, DiscretizeAddOnlyBnd) {
    auto A = test_bool_op(false, true);
    for (int i = 0; i < A.size(); ++i) {
        ASSERT_LT(A.type(i), 0) << format("Node %d (%s) is not a boundary node.", i, A.pos(i));
    }
}

TEST(Domains, DiscretizeAddFullDomain) {
    auto A = test_bool_op(true, true);
}


TEST(Domains, DiscretizeSubOnlyBnd) {
    auto A = test_bool_op(false, false);
    for (int i = 0; i < A.size(); ++i) {
        ASSERT_LT(A.type(i), 0) << format("Node %d (%s) is not a boundary node.", i, A.pos(i));
    }
}

TEST(Domains, DiscretizeSubFullDomain) {
    auto A = test_bool_op(true, false);
}


TEST(Domains, DISABLED_DomainDiscretizationUsageExample) {
    /// [Domain discretization usage example]
    BoxShape<Vec2d> b({0.0, 0.3}, {0.4, 2.7});
    DomainDiscretization<Vec2d> d = b.discretizeWithStep(0.1);
    std::cout << d.dim << std::endl;  // dimension of the domain
    d.size();  // number of points
    int i = 0;
    d.pos(i);  // position of node i
    if (d.type(i) < 0) {
        d.normal(i);  // boundary nodes have normals
    }
    d.boundary();  // indices of boundary nodes
    d.interior();  // indices of interior nodes
    d.shape();  // get geometric shape of the domain
    d.addBoundaryNode({0.001, 0.3}, -2, {-1, 0});  // add new boundary point
    d.addInternalNode({0.25, 0.75}, 2);  // add new internal point

    BallShape<Vec2d> c(0.0, 1.0);
    auto d2 = c.discretizeBoundaryWithStep(0.1);
    d -= d2;  // subtract or add discretizations

    std::cout << d << std::endl;
    /// [Domain discretization usage example]
}

TEST(Domains, DISABLED_DomainDiscretizationSubtractUsageExample) {
    /// [subtract example]
    double h = 0.1;
    DomainDiscretization<Vec2d> rect = BoxShape<Vec2d>({-1, -1}, {1, 1}).discretizeWithStep(h);
    DomainDiscretization<Vec2d> circ = BallShape<Vec2d>({0, 0}, 0.7).discretizeBoundaryWithStep(h);
    circ.types()[circ.boundary()] = -5;   // -1 to -4 taken by rectangle sides
    rect -= circ;
    std::cout << (rect.types() == -5).size() << std::endl;  // nonzero number is printed
    /// [subtract example]
}

TEST(Domains, DISABLED_DomainDiscretizationAddUsageExample) {
    /// [add example]
    double h = 0.1;
    DomainDiscretization<Vec2d> rect = BoxShape<Vec2d>({-1, -1}, {1, 1}).discretizeWithStep(h);
    rect.types()[rect.interior()] = 2;
    DomainDiscretization<Vec2d> circ = BallShape<Vec2d>({0, 0}, 0.7).discretizeWithStep(h);
    circ.types()[circ.boundary()] = -5;   // -1 to -4 taken by rectangle sides
    circ += rect;
    std::cout << (circ.types() == 2).size() << std::endl;  // nonzero number is printed
    /// [add example]
}

}  // namespace mm
