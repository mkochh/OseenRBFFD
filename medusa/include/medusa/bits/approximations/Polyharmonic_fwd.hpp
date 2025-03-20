#ifndef MEDUSA_BITS_APPROXIMATIONS_POLYHARMONIC_FWD_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_POLYHARMONIC_FWD_HPP_

/**
 * @file
 * Declaration of Polyharmonic RBF.
 *
 * @example test/approximations/Polyharmonic_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/utils/numutils.hpp>
#include <cmath>
#include "Operators_fwd.hpp"

namespace mm {

/**
 * Polyharmonic Radial Basis Function of odd order.
 * Specifically, this class represents a function
 * @f$ f(r^2) = r^k, k = 1, 3, \ldots@f$
 * as a function of squared distance with exponent @f$k@f$.
 * This class satisfies the @ref rbf-concept.
 *
 * @tparam k Order of the RBF. Must be odd.
 *
 * Usage example:
 * @snippet approximations/Polyharmonic_test.cpp Polyharmonic usage example
 * @ingroup approximations
 */
template <typename scal_t, int k = -1>
class Polyharmonic {
  public:
    static_assert((k > 0 && k % 2 == 1) || k == -1, "k must be odd or -1 (Dynamic)");
    typedef scal_t scalar_t;  ///< Scalar type used for computations.

  private:
    const int order_;  ///< Exponent of the RBF.

  public:
    /// Default constructor. Only applicable when order given as template argument.
    Polyharmonic();
    /// If template argument was -1, the order is accepted as runtime parameter.
    Polyharmonic(int order);

    /// Get order of the RBF.
    int order() const { return order_; }

    /**
     * Evaluate derivative of this RBF wrt. `r2` at given point.
     *
     * @param r2 Squared radial distance.
     * @param derivative Order of derivative to evaluate.
     * @throw Assertion fails if derivative order is negative.
     */
    inline scalar_t operator()(scalar_t r2, int derivative) const;

    /**
     * Evaluate Laplacian of this RBF wrt. `r` at given point.
     * @param r2 Squared radial distance.
     * @param lap Laplacian operator object.
     */
    template <int dimension>
    inline scalar_t operator()(scalar_t r2, Lap<dimension> lap) const;

    /// Evaluate this RBF given squared radial distance.
    inline scalar_t operator()(scalar_t r2) const { return ipow(std::sqrt(r2), order_); }

    /// Output basic information about given Gaussian RBF.
    template <typename S, int K>
    friend std::ostream& operator<<(std::ostream& os, const Polyharmonic<S, K>& m);
};

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_POLYHARMONIC_FWD_HPP_
