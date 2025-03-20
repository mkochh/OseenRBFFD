#ifndef MEDUSA_BITS_APPROXIMATIONS_POLYHARMONIC_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_POLYHARMONIC_HPP_

/**
 * @file
 * Implementation of Polyharmonic RBF.
 */

#include "Polyharmonic_fwd.hpp"
#include <medusa/bits/utils/numutils.hpp>

namespace mm {

template <typename scal_t, int k>
Polyharmonic<scal_t, k>::Polyharmonic() : order_(k) {
    assert_msg(order_ > 0, "Order must be supplied if compile type order is -1.");
}

template <typename scal_t, int k>
Polyharmonic<scal_t, k>::Polyharmonic(int order) : order_((k < 0) ? order : k) {
    assert_msg(order % 2 == 1, "Order must be odd, got %d.", order);
}

template <typename scal_t, int k>
scal_t Polyharmonic<scal_t, k>::operator()(scal_t r2, int derivative) const {
    assert_msg(derivative >= 0, "Derivative of negative order %d requested.", derivative);
    if (r2 < 1e-15) {
        return 0;
    }
    scalar_t r = std::sqrt(ipowneg(r2, order_ - 2*derivative));
    while (--derivative >= 0) {
        r *= (order_/2.0 - derivative);
    }
    return r;
}
/// @cond
template <typename scal_t, int k>
template <int dimension>
scal_t Polyharmonic<scal_t, k>::operator()(scal_t r2, Lap<dimension>) const {
    assert_msg(order_ > 2, "Order must be > 2 to compute the Laplacian, got %d.", k);
    return (dimension + order_ - 2)*order_*ipow(std::sqrt(r2), order_ - 2);
}
/// @endcond

/// Output basic information about given Gaussian RBF.
template <class S, int K>
std::ostream& operator<<(std::ostream& os, const Polyharmonic<S, K>& b) {
    return os << "Polyharmonic RBF of order " << b.order_ << '.';
}

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_POLYHARMONIC_HPP_
