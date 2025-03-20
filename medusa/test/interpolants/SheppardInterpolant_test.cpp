#include <medusa/bits/types/Vec.hpp>
#include <medusa/bits/utils/numutils.hpp>

#include "gtest/gtest.h"
#include "medusa/bits/interpolants/Sheppard_fwd.hpp"

namespace mm {

TEST(Interpolants, SheppardInterpolant1D) {
    /// [Sheppard interpolant usage example]
    int N = 4;                // Number of points.
    int num_closest = 2;      // Closest neighbors.
    Range<Vec1d> pos(N);      // Positions.
    Range<double> values(N);  // Values.

    for (int i = 0; i < N; i++) {
        pos[i] = i;
        values[i] = i;
    }

    // Sheppard's Interpolation.
    SheppardInterpolant<Vec1d, double> interpolant(pos, values);
    Vec1d point;
    point.setConstant(1.5);

    // Interpolated value.
    double value = interpolant(point, num_closest);
    /// [Sheppard interpolant usage example]

    EXPECT_NEAR(1.5, value, 1e-15);
}

TEST(Interpolants, SheppardInterpolant2D) {
    int N = 4;                // Number of points.
    Range<Vec2d> pos(N);      // Positions.
    Range<double> values(N);  // Values.

    // Rectangle corners.
    pos[0] = Vec2d{0, 0};
    pos[1] = Vec2d{1, 0};
    pos[2] = Vec2d{1, 1};
    pos[3] = Vec2d{0, 1};

    // Values in corners.
    for (int i = 0; i < N; i++) {
        values[i] = i;
    }

    // Sheppard's Interpolation.
    SheppardInterpolant<Vec2d, double> interpolant(pos, values);

    // Test points.
    int N_points = 30;
    auto x = Eigen::ArrayXf::LinSpaced(N_points, 0, 1);
    auto y = x;
    // Exclude corners.
    for (int i = 1; i < N_points - 1; i++) {
        for (int j = 1; j < N_points - 1; j++) {
            auto point = Vec2d{x(i), y(j)};
            // Interpolant value.
            double value = interpolant(point, N);
            // Expected value.
            Range<double> w(N);
            double w_sum = 0;
            for (int k = 0; k < N; k++) {
                w[k] = 1 / (point - pos[k]).squaredNorm();
                w_sum += w[k];
            }
            double expected_value = 0;
            for (int k = 0; k < N; k++) {
                expected_value += w[k] * values[k];
            }
            expected_value /= w_sum;

            EXPECT_NEAR(expected_value, value, 1e-15);
        }
    }
    // Test corner values.
    for (int i = 0; i < N; i++) {
        double value = interpolant(pos[i], N);

        EXPECT_NEAR(values[i], value, 1e-15);
    }
}

TEST(Interpolants, SheppardInterpolant3D) {
    int N = 3;                // Number of points.
    Range<Vec3d> pos(N);      // Positions.
    Range<double> values(N);  // Values.

    // Tetrahedron side corners.
    pos[0] = Vec3d{1, 0, 0};
    pos[1] = Vec3d{0, 1, 0};
    pos[2] = Vec3d{0, 0, 1};

    // Values in corners.
    for (int i = 0; i < N; i++) {
        values[i] = i;
    }

    // Sheppard's Interpolation.
    SheppardInterpolant<Vec3d, double> interpolant(pos, values);

    // Interpolated value at tetrahedron top.
    Vec3d point{0, 0, 0};
    double value = interpolant(point, 3);
    EXPECT_NEAR(1, value, 1e-15);

    // Interpolated value on the side.
    point = 0.2 * (pos[1] + pos[2]);
    value = interpolant(point, 2);
    EXPECT_NEAR(1.5, value, 1e-15);
}

TEST(Interpolants, SheppardInterpolantPower) {
    int N = 2;                // Number of points.
    int num_closest = 2;      // Closest neighbors.
    Range<Vec1d> pos(N);      // Positions.
    Range<double> values(N);  // Values.

    for (int i = 0; i < N; i++) {
        pos[i] = i;
        values[i] = i + 2;
    }

    // Sheppard's Interpolation.
    SheppardInterpolant<Vec1d, double> interpolant(pos, values);
    Vec1d point{0.75};

    // Interpolated value.
    for (int power = 1; power < 10; power++) {
        double value = interpolant(point, num_closest, power);

        // Expected value.
        Range<double> w(N);
        double w_sum = 0;
        for (int k = 0; k < N; k++) {
            w[k] = 1.0 / std::pow((point - pos[k]).norm(), power);
            w_sum += w[k];
        }
        double expected_value = 0;
        for (int k = 0; k < N; k++) {
            expected_value += w[k] * values[k];
        }
        expected_value /= w_sum;

        EXPECT_NEAR(expected_value, value, 1e-15);
    }
}

TEST(Interpolants, SheppardInterpolantSingleNeighbor) {
    int N = 2;                // Number of points.
    int num_closest = 1;      // Closest neighbors.
    Range<Vec1d> pos(N);      // Positions.
    Range<double> values(N);  // Values.

    for (int i = 0; i < N; i++) {
        pos[i] = 5 * i;
        values[i] = i + 2;
    }

    // Sheppard's Interpolation.
    SheppardInterpolant<Vec1d, double> interpolant(pos, values);

    // Interpolated value.
    Vec1d point{2.0};
    double value = interpolant(point, num_closest);

    EXPECT_NEAR(values[0], value, 1e-15);

    // Interpolated value.
    point.setConstant(2.6);
    value = interpolant(point, num_closest);

    EXPECT_NEAR(values[1], value, 1e-15);
}

TEST(Interpolants, SheppardInterpolantRegularization) {
    int N = 2;                // Number of points.
    int num_closest = 2;      // Closest neighbors.
    Range<Vec1d> pos(N);      // Positions.
    Range<double> values(N);  // Values.

    for (int i = 0; i < N; i++) {
        pos[i] = i;
        values[i] = i + 2;
    }

    // Sheppard's Interpolation.
    SheppardInterpolant<Vec1d, double> interpolant(pos, values);
    Vec1d point{0.75};

    // Interpolated value.
    for (double reg = 0.0; reg < 20.0; reg += 3.0) {
        double value = interpolant(point, num_closest, 2, reg);

        // Expected value.
        Range<double> w(N);
        double w_sum = 0;
        for (int k = 0; k < N; k++) {
            w[k] = 1.0 / ((point - pos[k]).squaredNorm() + reg);
            w_sum += w[k];
        }
        double expected_value = 0;
        for (int k = 0; k < N; k++) {
            expected_value += w[k] * values[k];
        }
        expected_value /= w_sum;

        EXPECT_NEAR(expected_value, value, 1e-15);
    }
}

TEST(Interpolants, SheppardInterpolantConfusionDistance) {
    int N = 2;                // Number of points.
    int num_closest = 2;      // Closest neighbors.
    Range<Vec1d> pos(N);      // Positions.
    Range<double> values(N);  // Values.

    for (int i = 0; i < N; i++) {
        pos[i] = i;
        values[i] = i + 2;
    }

    // Sheppard's Interpolation.
    SheppardInterpolant<Vec1d, double> interpolant(pos, values);

    // Interpolated value.
    for (double conf_dist = 0.1; conf_dist > 1e-14; conf_dist /= 10) {
        Vec1d point = Vec1d{conf_dist} - Vec1d{1e-15};
        double value = interpolant(point, num_closest, 2, 0, conf_dist);

        EXPECT_NEAR(values[0], value, 1e-15);
    }
}

}  // namespace mm
