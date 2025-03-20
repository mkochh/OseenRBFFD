#ifndef MEDUSA_BITS_APPROXIMATIONS_INVERSEMULTIQUADRIC_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_INVERSEMULTIQUADRIC_HPP_

/**
 * @file
 * Implementation of InverseMultiQuadratic RBF.
 */

#include "InverseMultiquadric_fwd.hpp"
#include <cmath>
#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/utils/numutils.hpp>

namespace mm {

template <typename scal_t>
InverseMultiquadric<scal_t>::InverseMultiquadric(scal_t shape) : shape_(shape) {
    assert_msg(shape_ > 0, "Shape should be greater than 0, got %s.", shape_);
}

template <class scal_t>
scal_t InverseMultiquadric<scal_t>::operator()(scal_t r2, int derivative) const {
    assert_msg(derivative >= 0, "Derivative of negative order %d requested.", derivative);
    scalar_t f = r2/shape_/shape_;
    scalar_t c = -0.5;
    for (int i = 1; i < derivative; ++i) {
        c *= (-0.5 - i);
    }
    return c / ipow(shape_, 2*derivative) / std::sqrt(ipow(f+1, 2*derivative + 1));
}

/// @cond
template <class scal_t>
template <int dimension>
scal_t InverseMultiquadric<scal_t>::operator()(scal_t r2, Lap<dimension>) const {
    scalar_t f = 1.0/shape_/shape_;
    scal_t inverse = 1.0 / std::sqrt(1+f*r2);
    return - dimension*f*ipow(inverse, 3) + 3*r2*f*f*ipow(inverse, 5);
}
/// @endcond

template <class scal_t>
scal_t InverseMultiquadric<scal_t>::operator()(scal_t r2) const {
    return 1.0 / std::sqrt(r2/shape_/shape_ + 1);
}

/// Output basic information about given basis function.
template <class S>
std::ostream& operator<<(std::ostream& os, const InverseMultiquadric<S>& b) {
    return os << "InverseMultiquadric RBFs with shape " << b.shape();
}

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_INVERSEMULTIQUADRIC_HPP_
