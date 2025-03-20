#ifndef MEDUSA_BITS_TYPES_VECTORFIELD_HPP_
#define MEDUSA_BITS_TYPES_VECTORFIELD_HPP_

/**
 * @file
 * Implementation of VectorField.
 */

#include "VectorField_fwd.hpp"

namespace mm {

template <typename Scalar, int dimension>
template <typename OtherDerived>
VectorField<Scalar, dimension>::VectorField(const Eigen::MatrixBase<OtherDerived>& other)
        : Base(other) {}

/// @cond
template <typename Scalar, int dimension>
template <typename OtherDerived>
VectorField<Scalar, dimension>& VectorField<Scalar, dimension>::operator=(
        const Eigen::MatrixBase<OtherDerived>& other) {
    Base::operator=(other); return *this;
}
/// @endcond

template <typename Scalar, int dimension>
VectorField<Scalar, dimension>&
VectorField<Scalar, dimension>::operator=(const Vec<Scalar, dimension>& v) {
    Base::rowwise() = v.transpose();
    return *this;
}

/// @cond
template <typename Scalar, int dimension>
template <typename OtherDerived>
VectorField<Scalar, dimension> VectorField<Scalar, dimension>::fromLinear(
        const Eigen::PlainObjectBase<OtherDerived>& other) {
    assert_msg(other.size() % dimension == 0, "The length of the given data is not a "
                                              "multiple of %d.", dimension);
    return Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, dimension>, 0>(
            other.data(), other.size()/dimension, dimension);
}
/// @endcond

}  // namespace mm

#endif  // MEDUSA_BITS_TYPES_VECTORFIELD_HPP_
