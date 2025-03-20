#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/domains/BallShape.hpp>
#include <medusa/bits/io/HDF.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(Domains, BoxShapeBegEndSwitch) {
    BoxShape<Vec1d> box1d({1.0}, {-1.0});
    BoxShape<Vec2d> box2d({1, -1}, {-1, 1});
    BoxShape<Vec3d> box3d({1, -1, 1}, {-1, 1, -1});

    for (int i = 0; i < box1d.dim; ++i) {
        EXPECT_DOUBLE_EQ(box1d.beg()[i], -1);
        EXPECT_DOUBLE_EQ(box1d.end()[i], 1);
    }
    for (int i = 0; i < box2d.dim; ++i) {
        EXPECT_DOUBLE_EQ(box2d.beg()[i], -1);
        EXPECT_DOUBLE_EQ(box2d.end()[i], 1);
    }
    for (int i = 0; i < box2d.dim; ++i) {
        EXPECT_DOUBLE_EQ(box2d.beg()[i], -1);
        EXPECT_DOUBLE_EQ(box2d.end()[i], 1);
    }
}

TEST(Domains, BoxShapeContainsBasic1D) {
    BoxShape<Vec1d> box1d({1.0}, {2.0});

    EXPECT_TRUE(box1d.contains({1.5}));
    EXPECT_TRUE(box1d.contains({1.01}));
    EXPECT_TRUE(box1d.contains({1.99}));
    EXPECT_TRUE(box1d.contains({1}));
    EXPECT_TRUE(box1d.contains({2}));
    EXPECT_TRUE(box1d.contains({1 - 1e-13}));
    EXPECT_TRUE(box1d.contains({2 + 1e-13}));
    EXPECT_FALSE(box1d.contains({3}));
    EXPECT_FALSE(box1d.contains({0}));
    EXPECT_FALSE(box1d.contains({-1}));
}

TEST(Domains, BoxShapeContainsBasic2D) {
    BoxShape<Vec2d> box2d({1, 3}, {1.2, 3.5});

    EXPECT_TRUE(box2d.contains({1.1, 3.4}));
    EXPECT_FALSE(box2d.contains({0.9, 3.4}));
    EXPECT_FALSE(box2d.contains({1.1, 2.3}));
    EXPECT_FALSE(box2d.contains({0.6, -1.4}));
}

TEST(Domains, BoxShapeContainsBasic3D) {
    BoxShape<Vec3d> box3d({0, 1, 2}, {-1, 2, 3});

    EXPECT_TRUE(box3d.contains({-0.5, 1.5, 2.5}));
    EXPECT_TRUE(box3d.contains({0, 1, 2}));
    EXPECT_FALSE(box3d.contains({-1.5, 1.5, 2.5}));
    EXPECT_FALSE(box3d.contains({-0.5, 2.5, 2.5}));
    EXPECT_FALSE(box3d.contains({-0.5, 1.5, 3.5}));
}

template<typename vec_t>
void testRectangleNormals(const DomainDiscretization<vec_t>& discretization,
                          const BoxShape<vec_t>& domain) {
    double tol = 1e-15;
    for (int i : discretization.types() < 0) {
        int on_bnd = 0;
        for (int j = 0; j < discretization.dim; ++j) {
            if (std::abs(discretization.positions()[i][j] - domain.beg()[j]) < tol) on_bnd++;
            if (std::abs(discretization.positions()[i][j] - domain.end()[j]) < tol) on_bnd++;
        }
        if (on_bnd == 1) {  // uniquely defined normal
            for (int j = 0; j < discretization.dim; ++j) {
                if (std::abs(discretization.positions()[i][j] - domain.beg()[j]) < tol) {
                    EXPECT_EQ(-1, discretization.normal(i)[j]);
                } else if (std::abs(discretization.positions()[i][j] - domain.end()[j]) < tol) {
                    EXPECT_EQ(1, discretization.normal(i)[j]);
                } else {
                    EXPECT_EQ(0, discretization.normal(i)[j]);
                }
            }
        }
    }
}

TEST(Domains, BoxShapeNormalsUniform) {
    BoxShape<Vec1d> domain1d({1}, {2});
    auto discretization1d = domain1d.discretizeBoundary({3});
    testRectangleNormals(discretization1d, domain1d);

    BoxShape<Vec2d> domain2d({1, 3}, {1.2, 3.5});
    auto discretization2d = domain2d.discretizeBoundary({3, 3});
    testRectangleNormals(discretization2d, domain2d);

    BoxShape<Vec3d> domain3d({0, 1, 2}, {-1, 2, 3});
    auto discretization3d = domain3d.discretizeBoundary({3, 3, 3});
    testRectangleNormals(discretization3d, domain3d);
}

TEST(Domains, BoxShapeleNormalsRandom) {
    BoxShape<Vec1d> domain1d({1}, {2});
    auto discretization1d = domain1d.discretizeWithStep(0.1);
    testRectangleNormals(discretization1d, domain1d);

    BoxShape<Vec2d> domain2d({1, 3}, {1.2, 3.5});
    auto discretization2d = domain2d.discretizeWithStep(0.1);
    testRectangleNormals(discretization2d, domain2d);

    BoxShape<Vec3d> domain3d({0, 1, 2}, {-1, 2, 3});
    auto discretization3d = domain3d.discretizeWithStep(0.1);
    testRectangleNormals(discretization3d, domain3d);
}

TEST(Domains, BoxDiscretize1d) {
    BoxShape<Vec1d> shape({-1.2}, {-0.4});
    auto d = shape.discretizeWithStep(0.1);
    Range<Vec1d> internal = d.positions()[d.interior()];
    std::sort(internal.begin(), internal.end());
    Range<double> expected = {-1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5};

    double tol = 1e-15;
    ASSERT_EQ(expected.size(), internal.size());
    int n = expected.size();
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(expected[i], internal[i][0], tol);
    }
}

TEST(Domains, BoxDiscretize2d) {
    BoxShape<Vec2d> shape({1, 3}, {1.4, 3.5});
    auto d = shape.discretizeWithStep(0.1);
    Range<Vec2d> internal = d.positions()[d.interior()];
    std::sort(internal.begin(), internal.end());
    Range<Vec2d> expected = {{1.1, 3.1}, {1.1, 3.2}, {1.1, 3.3}, {1.1, 3.4}, {1.2, 3.1},
                             {1.2, 3.2}, {1.2, 3.3}, {1.2, 3.4}, {1.3, 3.1}, {1.3, 3.2},
                             {1.3, 3.3}, {1.3, 3.4}};
    double tol = 1e-15;
    ASSERT_EQ(expected.size(), internal.size());
    int n = expected.size();
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(expected[i][0], internal[i][0], tol);
        EXPECT_NEAR(expected[i][1], internal[i][1], tol);
    }
}


TEST(Domains, BoxDiscretize3d) {
    BoxShape<Vec3d> shape({1, 3, -2.4}, {1.4, 3.5, -2.6+1e-15});
    auto d = shape.discretizeWithStep(0.1);
    Range<Vec3d> internal = d.positions()[d.interior()];
    std::sort(internal.begin(), internal.end());
    Range<Vec3d> expected = {{1.1, 3.1, -2.5}, {1.1, 3.2, -2.5}, {1.1, 3.3, -2.5}, {1.1, 3.4, -2.5},
                             {1.2, 3.1, -2.5}, {1.2, 3.2, -2.5}, {1.2, 3.3, -2.5}, {1.2, 3.4, -2.5},
                             {1.3, 3.1, -2.5}, {1.3, 3.2, -2.5}, {1.3, 3.3, -2.5},
                             {1.3, 3.4, -2.5}};
    double tol = 1e-15;
    ASSERT_EQ(expected.size(), internal.size());
    int n = expected.size();
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(expected[i][0], internal[i][0], tol);
        EXPECT_NEAR(expected[i][1], internal[i][1], tol);
        EXPECT_NEAR(expected[i][2], internal[i][2], tol);
    }
}

TEST(Domains, BoxDiscretizeBoundary3d) {
BoxShape<Vec3d> domain3d(0, 1);
    double dx = 0.05;
    double tol = 0.05;
    double s = 0;
    double c = 0;

    auto discretization3d = domain3d.discretizeBoundaryWithStep(dx);

    KDTree<Vec3d> tree(discretization3d.positions());
    for (int i = 0; i  < discretization3d.size(); ++i) {
        Range<double> distances2 = std::get<1>(tree.query(discretization3d.pos(i), 2));
        s += std::sqrt(distances2[1]);
        c += 1;
    }

    EXPECT_NEAR(dx, s/c, dx*tol);  // 5% error tolerance
}

TEST(Domains, BoxDiscretizeBoundaryWithDensity2d) {
    BoxShape<Vec2d> domain2d(0, 1);
    double dx1 = 0.05;
    double dx2 = 0.01;
    double tol = 0.05;
    double s1 = 0;
    double c1 = 0;
    double s2 = 0;
    double c2 = 0;

    auto fill = [=](const Vec2d& p) { return (p[0] < 0.5) ? dx1 : dx2; };
    auto discretization2d = domain2d.discretizeBoundaryWithDensity(fill);

    KDTree<Vec2d> tree(discretization2d.positions());
    for (int i = 0; i  < discretization2d.size(); ++i) {
        if (discretization2d.pos(i, 0) < 0.5 - tol) {
            Range<double> distances2 = std::get<1>(tree.query(discretization2d.pos(i), 2));
            s1 += std::sqrt(distances2[1]);
            c1 += 1;
        } else if (discretization2d.pos(i, 0) > 0.5 + tol) {
            Range<double> distances2 = std::get<1>(tree.query(discretization2d.pos(i), 2));
            s2 += std::sqrt(distances2[1]);
            c2 += 1;
        }
    }

    EXPECT_NEAR(dx1, s1/c1, dx1*tol);  // 5% error tolerance
    EXPECT_NEAR(dx2, s2/c2, dx2*tol);
}

TEST(Domains, BoxDiscretizeBoundaryWithDensity3d) {
    BoxShape<Vec3d> domain3d(0, 1);
    double dx1 = 0.05;
    double dx2 = 0.11;
    double tol = 0.05;
    double s1 = 0;
    double c1 = 0;
    double s2 = 0;
    double c2 = 0;

    auto fill = [=](const Vec3d& p) { return (p[0] < 0.5) ? dx1 : dx2; };
    auto discretization3d = domain3d.discretizeBoundaryWithDensity(fill);

    KDTree<Vec3d> tree(discretization3d.positions());
    for (int i = 0; i  < discretization3d.size(); ++i) {
        if (discretization3d.pos(i, 0) < 0.5 - tol) {
            Range<double> distances2 = std::get<1>(tree.query(discretization3d.pos(i), 2));
            s1 += std::sqrt(distances2[1]);
            c1 += 1;
        } else if (discretization3d.pos(i, 0) > 0.5 + tol) {
            Range<double> distances2 = std::get<1>(tree.query(discretization3d.pos(i), 2));
            s2 += std::sqrt(distances2[1]);
            c2 += 1;
        }
    }

    EXPECT_NEAR(dx1, s1/c1, dx1*tol);  // 5% error tolerance
    EXPECT_NEAR(dx2, s2/c2, dx2*tol);
}

TEST(Domains, DISABLED_BoxShapeUsageExample) {
    /// [BoxShape usage example]
    BoxShape<Vec3d> box({0, 0.3, -0.1}, {1, 2, 5});
    if (box.contains({0.5, 1.3, 3.4})) {
        // do something
    }
    auto d = box.discretizeWithStep(0.1);
    std::cout << box << std::endl;
    /// [BoxShape usage example]
    (void) d;  // remove unused warning
}

}  // namespace mm
