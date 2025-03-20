#include <medusa/bits/domains/discretization_helpers.hpp>
#include <medusa/bits/types/Vec_fwd.hpp>

#include "gtest/gtest.h"

namespace mm {


TEST(Domains, SphereDiscretization1d) {
    /// [SphereDiscretization usage example]
    std::mt19937 gen;
    double r = 2.3;
    const int dim = 1;
    auto v = discretization_helpers::SphereDiscretization<double, dim>::construct(r, 4, gen);
    /// [SphereDiscretization usage example]
    ASSERT_EQ(2, v.size());
    std::sort(v.begin(), v.end());
    EXPECT_EQ(-r, v[0][0]);
    EXPECT_EQ(r, v[1][0]);
}

TEST(Domains, SphereDiscretization2d) {
    double r = 2.3;
    auto v = discretization_helpers::SphereDiscretization<double, 2>::construct(r, 4);
    ASSERT_EQ(v.size(), 4);
    std::sort(v.begin(), v.end());
    double tol = 1e-15;
    EXPECT_NEAR(v[0][0], -r, tol);
    EXPECT_NEAR(v[0][1], 0, tol);
    EXPECT_NEAR(v[1][0], 0, tol);
    EXPECT_NEAR(v[1][1], -r, tol);
    EXPECT_NEAR(v[2][0], 0, tol);
    EXPECT_NEAR(v[2][1], r, tol);
    EXPECT_NEAR(v[3][0], r, tol);
    EXPECT_NEAR(v[3][1], 0, tol);
}

TEST(Domains, SphereDiscretization3d) {
    double r = 2.3;
    auto v = discretization_helpers::SphereDiscretization<double, 3>::construct(r, 6);
    decltype(v) expected = {{-1.15, -1.99186, 2.43932e-16}, {-1.15, -0.995929, -1.725},
                            {-1.15, -0.995929, 1.725}, {-1.15, 0.995929, -1.725},
                            {-1.15, 0.995929, 1.725}, {-1.15, 1.99186, 0},
                            {1.15, -1.99186, 2.43932e-16}, {1.15, -0.995929, -1.725},
                            {1.15, -0.995929, 1.725}, {1.15, 0.995929, -1.725},
                            {1.15, 0.995929, 1.725}, {1.15, 1.99186, 0}};
    std::sort(v.begin(), v.end());
    double tol = 1e-5;  // large due to low precision data above, but not really a problem
    ASSERT_EQ(expected.size(), v.size());
    int n = expected.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(expected[i][j], v[i][j], tol);
        }
    }
}

TEST(Domains, DiscretizeLineWithDensity) {
    /// [discretizeLineWithDensity]
    Vec2d p = 0.0, q = {1.0, 0.5};
    auto fn = [](const Vec2d& p) { return 0.1*(1+p[0]); };
    auto d = discretization_helpers::discretizeLineWithDensity(p, q, fn);
    /// [discretizeLineWithDensity]
    Range<Vec2d> expected = {{0.0894427, 0.0447214}, {0.186885, 0.0934427}, {0.293044, 0.146522},
                             {0.408697, 0.204349}, {0.534695, 0.267347}, {0.671962, 0.335981},
                             {0.821507, 0.410753}};
    // plot the points: endpoints are not included, spacing is variable
    double tol = 1e-6;
    ASSERT_EQ(expected.size(), d.size());
    int n = d.size();
    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(expected[i][0], d[i][0], tol);
        EXPECT_NEAR(expected[i][1], d[i][1], tol);
    }
}

}  // namespace mm
