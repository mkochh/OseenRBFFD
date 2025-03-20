#ifndef MEDUSA_BITS_TYPES_SCALARFIELD_HPP_
#define MEDUSA_BITS_TYPES_SCALARFIELD_HPP_

/**
 * @file
 * Implementation of ScalarField.
 */

#include "ScalarField_fwd.hpp"
#include <iostream>

namespace mm {

/// @cond
template <typename Scalar>
template <typename OtherDerived>
ScalarField<Scalar>::ScalarField(const Eigen::MatrixBase<OtherDerived>& other) : Base(other) {}

template <typename Scalar>
template <typename OtherDerived>
ScalarField<Scalar>& ScalarField<Scalar>::operator=(
        const Eigen::MatrixBase<OtherDerived>& other) {
    Base::operator=(other); return *this;
}
/// @endcond

}  // namespace mm

#endif  // MEDUSA_BITS_TYPES_SCALARFIELD_HPP_
