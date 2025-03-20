#include "medusa/bits/domains/cad_helpers.hpp"

/**
 * @file
 * Instantiations of commonly used cad helpers functions.
 */

/// @cond
template mm::Vec<double, 3> mm::cad_helpers::evaluate_b_spline(double t, int p,
        const Range<Vec<double, 3>>& control_points, const Range<double>& knots, int k);
template mm::Vec<double, 2> mm::cad_helpers::evaluate_b_spline(double t, int p,
        const Range<Vec<double, 2>>& control_points, const Range<double>& knots, int k);

template void mm::cad_helpers::generate_b_spline_derivative(int p,
        const Range<Vec<double, 2>>& control_points, const Range<double>& knots,
        Range<Vec<double, 2>>& der_control_points, Range<double>& der_knots);
template void mm::cad_helpers::generate_b_spline_derivative(int p,
        const Range<Vec<double, 3>>& control_points, const Range<double>& knots,
        Range<Vec<double, 3>>& der_control_points, Range<double>& der_knots);

template void mm::cad_helpers::generate_b_spline_derivative_control_points(int p,
        const Range<Vec<double, 2>>& control_points, const Range<double>& knots,
        Range<Vec<double, 2>>& der_control_points);
template void mm::cad_helpers::generate_b_spline_derivative_control_points(int p,
        const Range<Vec<double, 3>>& control_points, const Range<double>& knots,
        Range<Vec<double, 3>>& der_control_points);

template void mm::cad_helpers::generate_b_spline_derivative_knots(const Range<double>& knots,
                                                                  Range<double>& der_knots);
/// @endcond
