#include <medusa/bits/types/VectorField.hpp>
#include <medusa/bits/types/Vec.hpp>
#include <vector>

#include "gtest/gtest.h"

namespace mm {

TEST(Types, VectorFieldConstruct) {
    int n = 5;
    VectorField2d v(n);
    v = -3.4;
    for (int i = 0; i < n; ++i) {
        EXPECT_EQ(-3.4, v(i, 0));
        EXPECT_EQ(-3.4, v(i, 1));
    }

    VectorField2d a(v.rows());
    EXPECT_EQ(a.rows()*a.cols(), a.size());

    VectorField2d b(v);
    EXPECT_TRUE((b-v).isZero(0));

    std::vector<int> x;
    VectorField2d c(x.size());
    EXPECT_EQ(0, x.size());


    v = Vec2d(4.1, 5.5);
    for (int i = 0; i < n; ++i) {
        EXPECT_EQ(4.1, v(i, 0));
        EXPECT_EQ(5.5, v(i, 1));
    }

    // componentwise access
    v.c(1) = v.c(0);
    for (int i = 0; i < n; ++i) {
        EXPECT_EQ(4.1, v(i, 0));
        EXPECT_EQ(4.1, v(i, 1));
    }
}

TEST(Types, VectorFieldAccess) {
    int n = 5;
    VectorField2d v(n);
    v.setZero();

    v[2][1] = 4;
    for (int i = 0; i < n; ++i) {
        if (i == 2) {
            EXPECT_EQ(0.0, v(i, 0));
            EXPECT_EQ(4.0, v(i, 1));
        } else {
            EXPECT_EQ(0.0, v(i, 0));
            EXPECT_EQ(0.0, v(i, 1));
        }
    }

    v(4, 0) = -5.2;
    for (int i = 0; i < n; ++i) {
        if (i == 2) {
            EXPECT_EQ(0.0, v(i, 0));
            EXPECT_EQ(4.0, v(i, 1));
        } else if (i == 4) {
            EXPECT_EQ(-5.2, v(i, 0));
            EXPECT_EQ(0.0, v(i, 1));
        } else {
            EXPECT_EQ(0.0, v(i, 0));
            EXPECT_EQ(0.0, v(i, 1));
        }
    }

    v[0] = Vec2d(2.4, -1.3);
    for (int i = 0; i < n; ++i) {
        if (i == 0) {
            EXPECT_EQ(2.4, v(i, 0));
            EXPECT_EQ(-1.3, v(i, 1));
        } else if (i == 2) {
            EXPECT_EQ(0.0, v(i, 0));
            EXPECT_EQ(4.0, v(i, 1));
        } else if (i == 4) {
            EXPECT_EQ(-5.2, v(i, 0));
            EXPECT_EQ(0.0, v(i, 1));
        } else {
            EXPECT_EQ(0.0, v(i, 0));
            EXPECT_EQ(0.0, v(i, 1));
        }
    }
}

TEST(Types, VectorFieldMuliIndexAccess) {
    int n = 5;
    VectorField3d v(n);
    v.setZero();

    indexes_t idx = {1, 2, 3};
    v[idx] = 4;
    for (int i = 0; i < n; ++i) {
        if (i == 1 || i == 2 || i == 3) {
            EXPECT_EQ(4.0, v(i, 0));
            EXPECT_EQ(4.0, v(i, 1));
            EXPECT_EQ(4.0, v(i, 2));
        } else {
            EXPECT_EQ(0.0, v(i, 0));
            EXPECT_EQ(0.0, v(i, 1));
            EXPECT_EQ(0.0, v(i, 2));
        }
    }

    v[idx] = {2.3, -4.55, 3.212};
    for (int i = 0; i < n; ++i) {
        if (i == 1 || i == 2 || i == 3) {
            EXPECT_EQ(2.3, v(i, 0));
            EXPECT_EQ(-4.55, v(i, 1));
            EXPECT_EQ(3.212, v(i, 2));
        } else {
            EXPECT_EQ(0.0, v(i, 0));
            EXPECT_EQ(0.0, v(i, 1));
            EXPECT_EQ(0.0, v(i, 2));
        }
    }

    v(idx).setConstant(5);
    for (int i = 0; i < n; ++i) {
        if (i == 1 || i == 2 || i == 3) {
            EXPECT_EQ(5.0, v(i, 0));
            EXPECT_EQ(5.0, v(i, 1));
            EXPECT_EQ(5.0, v(i, 2));
        } else {
            EXPECT_EQ(0.0, v(i, 0));
            EXPECT_EQ(0.0, v(i, 1));
            EXPECT_EQ(0.0, v(i, 2));
        }
    }
}

TEST(Types, VectorFieldElement) {
    VectorField3d v(5);
    v.setZero();
    EXPECT_EQ(1, v[3].cols());
    EXPECT_EQ(3, v[3].rows());
    Vec3d c(3, 4.2, -3.4);
    Vec3d a = v[2] + c;
    EXPECT_EQ(c, a);
}

TEST(Types, VectorFieldLinearization) {
    int n = 8;
    VectorField3d v(n);
    v.setRandom();
    Eigen::VectorXd x = v.asLinear();
    VectorField3d v2 = VectorField3d::fromLinear(x);
    EXPECT_TRUE((v2-v).isZero(0));
}

TEST(Types, VectorFieldIterate) {
    int n = 3;
    VectorField2d v(n);
    v.setRandom();
    std::vector<double> coef;
    for (int i = 0; i < n; ++i) coef.push_back(v(i, 0));
    for (int i = 0; i < n; ++i) coef.push_back(v(i, 1));

    std::vector<double> coef2;
    for (double x : v) coef2.push_back(x);

    EXPECT_EQ(coef, coef2);
}

TEST(Types, DISABLED_VectorFieldUsageExample) {
    /// [Vector field usage example]
    int n = 7;
    VectorField2d v(n);
    int s = v.size();  // returns the number of elements
    assert(s == v.dim*n);
    assert(n == v.rows());
    v = 0.0;  // set to constant value
    v = Vec2d(3.8, 5.6);  // set to constant vector
    v(4, 0) = -1.2;  // set one coefficient
    v[3][1] = 4.5;  // set one coefficient, alternative syntax
    v[4] = Vec2d(3.4, 5.2);  // set one vector
    v[{1, 3, 4}] = Vec2d(-4.5, 2.5);  // set a subfield to a constant vector
    v[{0, 6, 2}] = 0.0;  // set a subfield to a constant scalar
    auto c = v.c(0);  // get all first components
    // you can do any matrix operators with the vector field:
    v *= -2.4;
    v += v;
    std::cout << v << std::endl;
    // convert to and from linear representation
    Eigen::VectorXd x = v.asLinear();
    VectorField3d v2 = VectorField3d::fromLinear(x);  // equal to v
    /// [Vector field usage example]
    (void) s;  // unused
    (void) c;
}

}  // namespace mm
