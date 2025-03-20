#ifndef MEDUSA_BITS_DOMAINS_COMPUTE_NORMAL_HPP_
#define MEDUSA_BITS_DOMAINS_COMPUTE_NORMAL_HPP_

/**
 * @file
 * Implementation of normal computation.
 */

#include "compute_normal_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/LU>

namespace mm {
namespace surface_fill_internal {

/// @cond
template <typename scalar_t, int dim_from, int dim_to>
Vec<scalar_t, dim_to> compute_normal(Eigen::Matrix<scalar_t, dim_to, dim_from> J) {
    static_assert(dim_from < dim_to, "At least one free dimension must be present.");
    Eigen::JacobiSVD<decltype(J)> svd(J, Eigen::ComputeFullU);
    assert_msg(svd.rank() == dim_from, "Jacobi matrix does not have full rank.");
    Vec<scalar_t, dim_to> normal = svd.matrixU().col(dim_from);
    // Find correct orientation
    Eigen::Matrix<double, dim_from + 1, dim_from + 1> M;
    M.template topLeftCorner<dim_from + 1, dim_from>() =
            J.template topLeftCorner<dim_from + 1, dim_from>();
    M.template rightCols<1>() = normal.template tail<dim_from + 1>();
    normal *= (M.determinant() > 0) ? -1.0 : 1.0;
    return normal;
}
/// @endcond

}  // namespace surface_fill_internal
}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_COMPUTE_NORMAL_HPP_
