#ifndef MEDUSA_BITS_DOMAINS_CAD_HELPERS_HPP_
#define MEDUSA_BITS_DOMAINS_CAD_HELPERS_HPP_

/**
 * @file
 * Implementation of functions evaluating basis splines (B-splines) and other useful functions for CAD.
 */

#include "cad_helpers_fwd.hpp"
#include "medusa/bits/utils/assert.hpp"
#include <algorithm>

namespace mm {

namespace cad_helpers {

template <typename scalar_t, int dim>
Vec<scalar_t, dim> evaluate_b_spline(scalar_t t, int p,
        const Range<Vec<scalar_t, dim>>& control_points, const Range<scalar_t>& knots, int k) {
    // Check for special cases.
    if (t == *(knots.end() - 1)) {
        return *(control_points.end() - 1);
    } else if (t == *knots.begin()) {
        return *control_points.begin();
    }

    assert_msg(knots.size() >= k + p + 1 && k - p >= 0,
            "Not enough knots or parameter not within range. Try padding the knot vector p times.");

    // De Boor's algorithm.
    Range<Vec<scalar_t, dim>> d(control_points.begin() + (k - p), control_points.begin() + (k + 1));

    for (int r = 1; r < p + 1; r++) {
        for (int j = p; j > r - 1; j--) {
            scalar_t alpha = (t - knots[j + k - p]) / (knots[j + 1 + k - r] - knots[j + k - p]);
            d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j];
        }
    }

    return d[p];
}

template <typename scalar_t, int dim>
Vec<scalar_t, dim> evaluate_b_spline(scalar_t t, int p,
        const Range<Vec<scalar_t, dim>>& control_points, const Range<scalar_t>& knots,
        scalar_t epsilon) {
    if (t > knots[knots.size() - 1] && abs(t - knots[knots.size() - 1]) < epsilon) {
        t = knots[knots.size() - 1];
    } else if (t < knots[0] && abs(t - knots[0]) < epsilon) {
        t = knots[0];
    }

    // Binary search for k.
    int k = std::upper_bound(knots.begin(), knots.end(), t) - knots.begin();
    k--;

    return evaluate_b_spline(t, p, control_points, knots, k);
}

template <typename scalar_t, int dim>
void generate_b_spline_derivative(int p, const Range<Vec<scalar_t, dim>>& control_points,
        const Range<scalar_t>& knots, Range<Vec<scalar_t, dim>>& der_control_points,
        Range<scalar_t>& der_knots) {
    generate_b_spline_derivative_control_points(p, control_points, knots, der_control_points);
    generate_b_spline_derivative_knots(knots, der_knots);
}

template <typename scalar_t, int dim>
void generate_b_spline_derivative_control_points(int p,
        const Range<Vec<scalar_t, dim>>& control_points, const Range<scalar_t>& knots,
        Range<Vec<scalar_t, dim>>& der_control_points) {
    // Derivative has one less control point.
    der_control_points.resize(control_points.size() - 1);

    // Calculate control points.
    for (int i = 0; i < control_points.size() - 1; i++) {
        scalar_t temp = ((scalar_t)p) / (knots[i + p + 1] - knots[i + 1]);
        for (int j =0; j < dim; j++) {
            der_control_points[i](j) = temp * (control_points[i + 1](j) - control_points[i](j));
        }
    }
}

template <typename scalar_t>
void generate_b_spline_derivative_knots(const Range<scalar_t>& knots,
        Range<scalar_t>& der_knots) {
    // Derivative's order is one less than original order.
    der_knots.resize(knots.size() - 2);

    // Calculate derivative knots.
    for (int i = 0; i < der_knots.size(); i++) {
        der_knots[i] = knots[i + 1];
    }
}

}  // namespace cad_helpers
}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_CAD_HELPERS_HPP_
