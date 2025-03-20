#ifndef MEDUSA_BITS_DOMAINS_COMPUTE_NORMAL_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_COMPUTE_NORMAL_FWD_HPP_

/**
 * @file
 * Compute normals from Jacobian matrix.
 *
 * @example test/domains/compute_normal_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/types/Vec.hpp>

namespace mm {
/// Internal implementations needed by GeneralSurfaceFill.
namespace surface_fill_internal {

/**
 * Compute normal from given Jacobian.
 * @tparam scalar_t Numeric scalar type, e.g. `double` or `float`.
 * @tparam dim_from Parametric domain dimension.
 * @tparam dim_to Target domain dimension.
 * @param J `dim_to x dim_from` Jacobi matrix, which must have full rank.
 * @return A unit vector `n`, such that it is orthogonal to all columns of `J` and
 * that the determinant of `[J n]` is negative. This corresponds to the right-hand side
 * 90 deg rotation in 2d and the negative cross product in 3d.
 * If `dim_to == dim_from+1` this uniquely defines the normal.
 * @ingroup domains
 * @anchor compute_normal
 */
template <typename scalar_t, int dim_from, int dim_to>
Vec<scalar_t, dim_to> compute_normal(Eigen::Matrix<scalar_t, dim_to, dim_from> J);

}  // namespace surface_fill_internal
}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_COMPUTE_NORMAL_FWD_HPP_
