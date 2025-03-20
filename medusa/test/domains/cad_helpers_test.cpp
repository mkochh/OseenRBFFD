#include "medusa/bits/domains/cad_helpers.hpp"
#include "medusa/bits/types/Vec.hpp"
#include "medusa/bits/types/Range_fwd.hpp"
#include <cmath>

#include "gtest/gtest.h"

namespace mm {
namespace cad_helpers {

TEST(DomainEngines, B_splineEvaluate) {
    int p = 2;
    Range<double> knots({0, 0, 0, 1, 2, 3, 3, 4, 4, 4});
    Range<Vec2d> cp({Vec2d(0, 1), Vec2d(1, 1), Vec2d(3, 4),
                     Vec2d(4, 2), Vec2d(5, 3), Vec2d(6, 4), Vec2d(7, 3)});

    Vec2d pt = evaluate_b_spline(2.4, p, cp, knots, 4);
    EXPECT_LT((pt - Vec2d(3.98, 2.52)).norm(), 1e-15);

    pt = evaluate_b_spline(1.2, p, cp, knots, 3);
    EXPECT_LT((pt - Vec2d(2.38, 3)).norm(), 1e-15);

    pt = evaluate_b_spline(0.4, p, cp, knots);
    EXPECT_LT((pt - Vec2d(0.8, 1.24)).norm(), 1e-15);

    pt = evaluate_b_spline(3.9, p, cp, knots);
    EXPECT_LT((pt - Vec2d(6.8, 3.18)).norm(), 1e-15);

    pt = evaluate_b_spline(0.0, p, cp, knots);
    EXPECT_LT((pt - Vec2d(0, 1)).norm(), 1e-15);

    pt = evaluate_b_spline(4.0, p, cp, knots);
    EXPECT_LT((pt - Vec2d(7, 3)).norm(), 1e-15);
}

TEST(DomainEngines, B_splineDerivative) {
    int p = 2;
    Range<double> knots({0, 0, 0, 1, 2, 3, 3, 4, 4, 4});
    Range<Vec2d> cp({Vec2d(0, 1), Vec2d(1, 1), Vec2d(3, 4),
                     Vec2d(4, 2), Vec2d(5, 3), Vec2d(6, 4), Vec2d(7, 3)});

    Range<double> d_knots;
    Range<Vec2d> d_cp;
    cad_helpers::generate_b_spline_derivative(p, cp, knots, d_cp, d_knots);

    Vec2d pt = cad_helpers::evaluate_b_spline(3.5, p - 1, d_cp, d_knots);
    EXPECT_LT((pt - Vec2d(2, 0)).norm(), 1e-15);

    pt = cad_helpers::evaluate_b_spline(0.0, p - 1, d_cp, d_knots);
    EXPECT_LT((pt - Vec2d(2, 0)).norm(), 1e-15);

    pt = cad_helpers::evaluate_b_spline(1.5, p - 1, d_cp, d_knots);
    EXPECT_LT((pt - Vec2d(1.5, 0.5)).norm(), 1e-15);

    pt = cad_helpers::evaluate_b_spline(4.0, p - 1, d_cp, d_knots);
    EXPECT_LT((pt - Vec2d(2, -2)).norm(), 1e-15);

    pt = cad_helpers::evaluate_b_spline(2.5, p - 1, d_cp, d_knots);
    EXPECT_LT((pt - Vec2d(1.5, 0)).norm(), 1e-15);
}

}  // namespace cad_helpers
}  // namespace mm
