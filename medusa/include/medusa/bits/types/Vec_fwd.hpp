#ifndef MEDUSA_BITS_TYPES_VEC_FWD_HPP_
#define MEDUSA_BITS_TYPES_VEC_FWD_HPP_

/**
 * @file
 * Declaration of Vec.
 *
 * @example test/types/Vec_test.cpp
 * @sa MatrixAddons.hpp
 * @sa MatrixBaseAddons.hpp
 */

#include <medusa/Config.hpp>
#include <Eigen/Core>

namespace mm {

/**
 * Fixed size vector type, representing a mathematical 1d/2d/3d vector.
 * This is a special case of fixed size Eigen type. See
 * our `Matrix` and `MatrixBase` add-ons for more functionality.
 *
 * @sa MatrixAddons.hpp
 * @sa MatrixBaseAddons.hpp
 *
 * Usage example:
 * @snippet types/Vec_test.cpp Vec usage example
 * @ingroup types
 */
template <class scalar_t, int dim>
using Vec = Eigen::Matrix<scalar_t, dim, 1, Eigen::ColMajor|Eigen::AutoAlign, dim, 1>;

typedef Vec<double, 1> Vec1d;  ///< Convenience typedef for 1d vector of doubles.
typedef Vec<double, 2> Vec2d;  ///< Convenience typedef for 2d vector of doubles.
typedef Vec<double, 3> Vec3d;  ///< Convenience typedef for 3d vector of doubles.

typedef Vec<std::complex<double >, 1> Vec1cd;
///< Convenience typedef for 1d vector of complex doubles.
typedef Vec<std::complex<double >, 2> Vec2cd;
///< Convenience typedef for 2d vector of complex doubles.
typedef Vec<std::complex<double >, 3> Vec3cd;
///< Convenience typedef for 3d vector of complex doubles.
}  // namespace mm

#endif  // MEDUSA_BITS_TYPES_VEC_FWD_HPP_
