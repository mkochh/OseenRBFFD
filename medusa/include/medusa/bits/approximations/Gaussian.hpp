#ifndef MEDUSA_BITS_APPROXIMATIONS_GAUSSIAN_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_GAUSSIAN_HPP_

/**
 * @file
 * Implementation of Gaussian RBF.
 */

#include "Gaussian_fwd.hpp"
#include <cmath>
#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/utils/numutils.hpp>

namespace mm {

template <typename scal_t>
Gaussian<scal_t>::Gaussian(scal_t shape) : shape_(shape) {
    assert_msg(shape_ > 0, "Shape should be greater than 0, got %s.", shape_);
}

template <class scal_t>
scal_t Gaussian<scal_t>::operator()(scal_t r2, int derivative) const {
    assert_msg(derivative >= 0, "Derivative of negative order %d requested.", derivative);
    scal_t f = -1.0/shape_/shape_;
    return std::exp(r2*f) * ipow(f, derivative);
}
/// @cond
template <class scal_t>
template <int dimension>
scal_t Gaussian<scal_t>::operator()(scalar_t r2, Lap<dimension>) const {
    scal_t f = -1.0/shape_/shape_;
    return (2*dimension*f + 4*r2*(f*f)) * std::exp(r2*f);
}
/// @endcond
/// Output basic information about given Gaussian RBF.
template <class S>
std::ostream& operator<<(std::ostream& os, const Gaussian<S>& b) {
    return os << "Gaussian RBF with shape " << b.shape();
}

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_GAUSSIAN_HPP_
