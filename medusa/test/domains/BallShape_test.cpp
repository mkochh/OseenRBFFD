#include <medusa/bits/domains/BallShape.hpp>
#include <medusa/bits/domains/DomainShape.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(Domains, BallContainsBasic1D) {
    BallShape<Vec1d> cs1({3}, 1);

    EXPECT_TRUE(cs1.contains({3.5}));
    EXPECT_TRUE(cs1.contains({4}));
    EXPECT_FALSE(cs1.contains({1}));
    EXPECT_FALSE(cs1.contains({5}));
}

TEST(Domains, BallContainsBasic2D) {
    BallShape<Vec2d> cs2({-1, 2}, 5);

    EXPECT_TRUE(cs2.contains({-0.9, 1.8}));
    EXPECT_TRUE(cs2.contains({-1.1, 2.2}));
    EXPECT_TRUE(cs2.contains({0, 3}));
    EXPECT_TRUE(cs2.contains({4, 2}));
    EXPECT_TRUE(cs2.contains({4 + 1e-13, 2}));
    EXPECT_FALSE(cs2.contains({-8, 3}));
    EXPECT_FALSE(cs2.contains({12, 3}));
}

TEST(Domains, BallContainsBasic3D) {
    BallShape<Vec3d> cs3({4, -2, 0}, 4);

    EXPECT_TRUE(cs3.contains({4, 2, 0}));
    EXPECT_TRUE(cs3.contains({4, -2, 3}));
    EXPECT_FALSE(cs3.contains({10, 3, 4}));
    EXPECT_FALSE(cs3.contains({-5, 3, 0}));
}

TEST(Domains, BallDiscretizeBoundary1d) {
    BallShape<Vec1d> b(-1, 0.4);
    DomainDiscretization<Vec1d> d = b.discretizeBoundaryWithStep(0.1);
    auto p = d.positions();
    std::sort(p.begin(), p.end());
    EXPECT_EQ(-1.4, p[0][0]);
    EXPECT_EQ(-0.6, p[1][0]);
}

TEST(Domains, BallDiscretizeBoundary2d) {
    BallShape<Vec2d> c({-1, 2.2}, 1.7);
    auto d = c.discretizeBoundaryWithStep(0.5);
    auto p = d.positions();
    std::vector<Vec2d> expected = {
            {0.7, 2.2}, {0.6311380551446455, 2.678945346630431},
            {0.4301310058130081, 3.119089389674516}, {0.1132632477069846, 3.484774276402239},
            {-0.293794477896793, 3.746374392102681}, {-0.7580647749354152, 3.882696451197586},
            {-1.241935225064585, 3.882696451197586}, {-1.706205522103207, 3.746374392102681},
            {-2.113263247706985, 3.484774276402239}, {-2.430131005813008, 3.119089389674516},
            {-2.631138055144645, 2.678945346630431}, {-2.7, 2.2},
            {-2.631138055144646, 1.72105465336957}, {-2.430131005813008, 1.280910610325484},
            {-2.113263247706985, 0.9152257235977614}, {-1.706205522103207, 0.6536256078973188},
            {-1.241935225064585, 0.5173035488024147}, {-0.7580647749354157, 0.5173035488024145},
            {-0.2937944778967937, 0.6536256078973186}, {0.113263247706985, 0.9152257235977614},
            {0.4301310058130081, 1.280910610325485}, {0.6311380551446455, 1.72105465336957}};
    double tol = 1e-15;
    ASSERT_EQ(expected.size(), p.size());
    int n = expected.size();
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(expected[i][0], p[i][0], tol);
        EXPECT_NEAR(expected[i][1], p[i][1], tol);
    }
}

TEST(Domains, BallDiscretizeBoundary3d) {
    BallShape<Vec3d> domain3d(0, 1);
    double dx = 0.05;
    double tol = 0.05;
    double s = 0;
    double c = 0;

    auto discretization3d = domain3d.discretizeBoundaryWithStep(dx);
    KDTree<Vec3d> tree(discretization3d.positions());
    for (int i = 0; i  < discretization3d.size(); ++i) {
        Range<double> distances2 = tree.query(discretization3d.pos(i), 2).second;
        s += std::sqrt(distances2[1]);
        c += 1;
    }

    EXPECT_NEAR(dx, s/c, dx*tol);  // 5% error tolerance
}

TEST(Domains, BallDiscretizeBoundaryWithDensity2d) {
    BallShape<Vec2d> domain2d(0, 1);
    double dx1 = 0.05;
    double dx2 = 0.01;
    double tol = 0.05;
    double s1 = 0;
    double c1 = 0;
    double s2 = 0;
    double c2 = 0;

    auto fill = [=](const Vec2d& p) { return (p[0] < 0.0) ? dx1 : dx2; };
    auto discretization2d = domain2d.discretizeBoundaryWithDensity(fill);

    KDTree<Vec2d> tree(discretization2d.positions());
    for (int i = 0; i  < discretization2d.size(); ++i) {
        if (discretization2d.pos(i, 0) < 0.0 - tol) {
            Range<double> distances2 = tree.query(discretization2d.pos(i), 2).second;
            s1 += std::sqrt(distances2[1]);
            c1 += 1;
        } else if (discretization2d.pos(i, 0) > 0.0 + tol) {
            Range<double> distances2 = tree.query(discretization2d.pos(i), 2).second;
            s2 += std::sqrt(distances2[1]);
            c2 += 1;
        }
    }

    EXPECT_NEAR(dx1, s1/c1, dx1*tol);  // 5% error tolerance
    EXPECT_NEAR(dx2, s2/c2, dx2*tol);
}

TEST(Domains, BallDiscretizeBoundaryWithDensity3d) {
    // TODO(jureslak): this test fails if tolerance is 0.05, but it shouldn't.
    // Some generated points are too close.
    double tol = 0.1;
    double s1 = 0;
    double c1 = 0;
    double s2 = 0;
    double c2 = 0;
    double dx1 = 0.2;
    double dx2 = 0.1;

    BallShape<Vec3d> domain3d(0, 1);

    auto fill = [=](const Vec3d& p) { return (p[2] < 0.0) ? dx1 : dx2; };  // z-coordinate
    auto discretization3d = domain3d.discretizeBoundaryWithDensity(fill);

    // Make kd tree, iterate over the points on the boundary and average the distance between them.
    KDTree<Vec3d> tree(discretization3d.positions());
    for (int i = 0; i  < discretization3d.size(); ++i) {
        if (discretization3d.pos(i, 2) < 0.0 - tol) {
            Range<double> distances2 = tree.query(discretization3d.pos(i), 2).second;
            s1 += std::sqrt(distances2[1]);
            c1 += 1;
        } else if (discretization3d.pos(i, 2) > 0.0 + tol) {
            Range<double> distances2 = tree.query(discretization3d.pos(i), 2).second;
            s2 += std::sqrt(distances2[1]);
            c2 += 1;
        }
    }

    EXPECT_NEAR(dx1, s1/c1, dx1*tol);  // 5% error tolerance
    EXPECT_NEAR(dx2, s2/c2, dx2*tol);
}

TEST(Domains, BallShapeNormalsRandom1d) {
    BallShape<Vec1d> cd1({3}, 1);
    auto discr1d = cd1.discretizeWithStep(0.1);
    for (int i : discr1d.boundary()) {
        Vec1d n = (discr1d.pos(i) - cd1.center()).normalized();
        EXPECT_DOUBLE_EQ(n[0], discr1d.normal(i)[0]);
    }
}
TEST(Domains, BallShapeNormalsRandom2d) {
    BallShape<Vec2d> cd2({-1, 2}, 5);
    auto discr2d = cd2.discretizeWithStep(0.1);
    double tol = 1e-15;
    for (int i : discr2d.boundary()) {
        Vec2d n = (discr2d.pos(i) - cd2.center()).normalized();
        EXPECT_NEAR(n[0], discr2d.normal(i)[0], tol);
        EXPECT_NEAR(n[1], discr2d.normal(i)[1], tol);
    }
}
TEST(Domains, BallShapeNormalsRandom3d) {
    BallShape<Vec3d> cd3({4, -2, 0}, 4);
    auto discr3d = cd3.discretizeBoundaryWithStep(0.2);
    double tol = 1e-15;
    for (int i : discr3d.boundary()) {
        Vec3d n = (discr3d.pos(i) - cd3.center()).normalized();
        EXPECT_NEAR(n[0], discr3d.normal(i)[0], tol);
        EXPECT_NEAR(n[1], discr3d.normal(i)[1], tol);
        EXPECT_NEAR(n[2], discr3d.normal(i)[2], tol);
    }
}

TEST(Domains, BallDiscretize1d) {
    BallShape<Vec1d> b(-1, 0.4);
    DomainDiscretization<Vec1d> d = b.discretizeWithStep(0.1);
    auto p = d.positions();
    std::sort(p.begin(), p.end());
    std::vector<double> expected = {-1.4, -1.3, -1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6};
    double tol = 1e-15;
    ASSERT_EQ(expected.size(), p.size());
    int n = expected.size();
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(expected[i], p[i][0], tol);
    }
}

TEST(Domains, BallDiscretize2d) {
    BallShape<Vec2d> b({-1, 2.7}, 0.5);
    DomainDiscretization<Vec2d> d = b.discretizeWithStep(0.1);
    auto p = d.positions();
    int n = p.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            ASSERT_GE((p[i]-p[j]).norm(), 0.095);  // due to integer number of nodes
        }
    }
}

TEST(Domains, BallDiscretize4d) {
    BallShape<Vec<double, 4>> b(0, 0.5);
    EXPECT_DEATH(b.discretizeWithStep(0.1), "This domain does not support filling with density");

    auto fill = [=](const Vec<double, 4>&) { return 0.1; };
    EXPECT_DEATH(b.discretizeBoundaryWithDensity(fill),
                 "This domain does not support filling with density");
}

TEST(Domains, BallDscretizeBoundaryWithStep4d) {
    typedef Vec<double, 4> vec;
    BallShape<vec> b(0, 1);
    double dx = 0.15;  // Chosen to generate approximately 6000 points
    double tol = 0.05;
    double s = 0;
    double c = 0;

    auto d = b.discretizeBoundaryWithStep(dx);
    KDTree<vec> tree(d.positions());

    for (int i = 0; i  < d.size(); ++i) {
        Range<double> distances2 = tree.query(d.pos(i), 2).second;
        s += std::sqrt(distances2[1]);
        c += 1;
    }

    EXPECT_NEAR(dx, s/c, dx*tol);  // 5 % error tolerance
}

TEST(Domains, DISABLED_BallShapeUsageExample) {
    /// [BallShape usage example]
    BallShape<Vec2d> ball({-1, 0}, 7);
    if (ball.contains({0.5, 1.3})) {
        // do something
    }
    auto d = ball.discretizeWithStep(0.1);
    std::cout << ball << std::endl;
    /// [BallShape usage example]
    (void) d;  // remove unused warning
}

}  // namespace mm
