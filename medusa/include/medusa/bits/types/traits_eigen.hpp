#ifndef MEDUSA_BITS_TYPES_TRAITS_EIGEN_HPP_
#define MEDUSA_BITS_TYPES_TRAITS_EIGEN_HPP_

/**
 * @file
 * Type traits for vector and scalar fields for Eigen types.
 */

#include <medusa/Config.hpp>
#include "traits.hpp"
#include <Eigen/Core>

namespace mm {

/// The scalar_type trait definition for Eigen Matrix.
template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
struct scalar_type<Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>> {
    typedef _Scalar type;  ///< Underlying scalar type.
};

/// The vector_type trait definition for Eigen Matrix.
template <typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
struct vector_type<Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>> {
    typedef Eigen::Matrix<_Scalar, _Cols, 1> type;  ///< Underlying vector type.
};

}  // namespace mm

#endif  // MEDUSA_BITS_TYPES_TRAITS_EIGEN_HPP_
