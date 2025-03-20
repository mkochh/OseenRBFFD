#ifndef MEDUSA_BITS_DOMAINS_CAD_HELPERS_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_CAD_HELPERS_FWD_HPP_

/**
 * @file
 * Functions evaluating basis splines (B-splines) and other useful functions for CAD.
 *
 * @example test/domains/cad_helpers_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/types/Vec.hpp>
#include <medusa/bits/types/Range.hpp>

namespace mm {

/// Namespace for helper functions used in CAD.
namespace cad_helpers {

/**
 * Evaluate B-spline in one point using De Boor's algorithm - @f$\mathcal{O}(d p^2)@f$,
 * where @f$d@f$ is the number of dimensions.
 * Adds additional time complexity of @f$\mathcal{O}(\log |\mathrm{knots}|)@f$ if `k` is not given.
 * @tparam scalar_t Numeric scalar type, e.g. `double` or `float`.
 * @tparam dim Space dimension.
 * @param p B-spline order (starting from order 0).
 * @param t Point of evaluation.
 * @param control_points Range of control points.
 * @param knots Knot vector.
 * @param k Position of `t` in the knot vector, such that `control_points[k] <= t < control_points[k]`
 * @return Vector to the point on the B-spline corresponding to parameter `t`.
 */
template <typename scalar_t, int dim>
Vec<scalar_t, dim> evaluate_b_spline(scalar_t t, int p,
        const Range<Vec<scalar_t, dim>>& control_points, const Range<scalar_t>& knots, int k);

/**
 * Overload with binary search for position of `t` in `knots`.
 * If parameter is outside the domain by epsilon, translate it back to the edge of
 * the domain.
 */
template <typename scalar_t, int dim>
Vec<scalar_t, dim> evaluate_b_spline(scalar_t t, int p,
        const Range<Vec<scalar_t, dim>>& control_points, const Range<scalar_t>& knots,
        scalar_t epsilon = 1e-10);

/**
 * Generate control points and knot vector of a B-spline that is the first derivative
 * of the inputed B-spline. Note that its degree is `p - 1`.
 * @tparam scalar_t Numeric scalar type, e.g. `double` or `float`.
 * @tparam dim Space dimension.
 * @param p B-spline order (starting from order 0).
 * @param control_points Range of control points.
 * @param knots Knot vector.
 * @param der_control_points Range for derivative control points.
 * @param der_knots Range for derivative knots.
 */
template <typename scalar_t, int dim>
void generate_b_spline_derivative(int p, const Range<Vec<scalar_t, dim>>& control_points,
        const Range<scalar_t>& knots, Range<Vec<scalar_t, dim>>& der_control_points,
        Range<scalar_t>& der_knots);

/**
 * Generate control points of a B-spline that is the first derivative
 * of the inputed B-spline. Note that its degree is `p - 1`.
 * @tparam scalar_t Numeric scalar type, e.g. `double` or `float`.
 * @tparam dim Space dimension.
 * @param p B-spline order (starting from order 0).
 * @param control_points Range of control points.
 * @param knots Knot vector.
 * @param der_control_points Range for derivative control points.
 */
template <typename scalar_t, int dim>
void generate_b_spline_derivative_control_points(int p,
        const Range<Vec<scalar_t, dim>>& control_points, const Range<scalar_t>& knots,
        Range<Vec<scalar_t, dim>>& der_control_points);

/**
 * Generate knots of a B-spline that is the first derivative
 * of the inputed B-spline. Note that its degree is `p - 1`.
 * @tparam scalar_t Numeric scalar type, e.g. `double` or `float`.
 * @param knots Knot vector.
 * @param der_knots Range for derivative knots.
 */
template <typename scalar_t>
void generate_b_spline_derivative_knots(const Range<scalar_t>& knots,
        Range<scalar_t>& der_knots);

}  // namespace cad_helpers
}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_CAD_HELPERS_FWD_HPP_
