#include "medusa/bits/domains/compute_normal.hpp"

/**
 * @file
 * Instantiations of commonly used compute_normals instances.
 */

/// @cond
template mm::Vec<double, 3> mm::surface_fill_internal::compute_normal(Eigen::Matrix<double, 3, 2>);
template mm::Vec<double, 4> mm::surface_fill_internal::compute_normal(Eigen::Matrix<double, 4, 3>);
template mm::Vec<double, 3> mm::surface_fill_internal::compute_normal(Eigen::Matrix<double, 3, 1>);
template mm::Vec<double, 2> mm::surface_fill_internal::compute_normal(Eigen::Matrix<double, 2, 1>);
/// @endcond
