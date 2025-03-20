#ifndef MEDUSA_BITS_APPROXIMATIONS_WEIGHTFUNCTION_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_WEIGHTFUNCTION_HPP_

/**
 * @file
 * Implementations of weight functions.
 */

#include <cmath>
#include <medusa/bits/utils/assert.hpp>
#include "WeightFunction_fwd.hpp"

namespace mm {

/// Output info about given weight function.
template <typename V, typename R>
std::ostream& operator<<(std::ostream& os, const RBFWeight<V, R>& w) {
    return os << "RBFWeight constructed from " << w.rbf_;
}

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_WEIGHTFUNCTION_HPP_
