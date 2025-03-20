#ifndef MEDUSA_BITS_APPROXIMATIONS_INVERSEMULTIQUADRIC_FWD_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_INVERSEMULTIQUADRIC_FWD_HPP_

/**
 * @file
 * Declaration of InverseMultiQuadratic RBF.
 *
 * @example test/approximations/InverseMultiquadric_test.cpp
 */

#include <medusa/Config.hpp>
#include <iosfwd>
#include "Operators_fwd.hpp"

namespace mm {

/**
 * Inverse Multiquadric Radial Basis Function.
 * Specifically, this class represents a function
 * @f$\frac{1}{\sqrt{1 + r^2/\sigma^2}}@f$.
 * as a function of squared distance with shape parameter @f$\sigma @f$.
 * This class satisfies the @ref rbf-concept.
 * 
 * Usage example:
 * @snippet approximations/InverseMultiquadric_test.cpp InverseMultiquadric usage example
 * @ingroup approximations
 */
template <class scal_t>
class InverseMultiquadric {
  public:
    typedef scal_t scalar_t;  ///< Scalar type used for computations.

  private:
    scalar_t shape_;  ///< Shape parameter

  public:
    /// Creates a InverseMultiQuadratic RBF basis with shape parameter `shape`.
    InverseMultiquadric(scalar_t shape = 1);

    /// Returns shape parameter.
    scalar_t shape() const { return shape_ ;}

    /// Sets shape parameter to a new value.
    void setShape(scalar_t shape) {shape_ = shape;}

    /**
     * Evaluate derivative of this RBF wrt. `r2` at given point.
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

    /// Evaluate at a given point.
    inline scalar_t operator()(scalar_t r2) const;

    /// Output basic information about given basis function.
    template <typename V>
    friend std::ostream& operator<<(std::ostream& os, const InverseMultiquadric<V>& m);
};


}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_INVERSEMULTIQUADRIC_FWD_HPP_
