#ifndef MEDUSA_BITS_APPROXIMATIONS_GAUSSIAN_FWD_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_GAUSSIAN_FWD_HPP_

/**
 * @file
 * Declaration of Gaussian RBF.
 *
 * @example test/approximations/Gaussian_test.cpp
 */

#include <medusa/Config.hpp>
#include <iosfwd>
#include <cmath>
#include "Operators_fwd.hpp"
namespace mm {

/**
 * Gaussian Radial Basis Function.
 * Specifically, this class represents a function
 * @f$ f(r^2) = \exp(-r^2/\sigma^2) @f$
 * as a function of squared distance with shape parameter @f$\sigma @f$.
 * This class satisfies the @ref rbf-concept.
 *
 * Usage example:
 * @snippet approximations/Gaussian_test.cpp Gaussian usage example
 * @ingroup approximations
 */
template <class scal_t>
class Gaussian {
  public:
    typedef scal_t scalar_t;  ///< Scalar type used for computations.

  private:
    scalar_t shape_;  ///< Shape parameter.

  public:
    /// Creates a Gaussian RBF with shape parameter `shape`. The shape should be positive.
    Gaussian(scalar_t shape = 1.0);

    /// Returns shape parameter.
    scalar_t shape() const { return shape_; }

    /// Sets shape parameter to a new value.
    void setShape(scalar_t shape) { shape_ = shape; }

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

    /// Evaluate this RBF given squared radial distance.
    inline scalar_t operator()(scalar_t r2) const { return std::exp(-r2/shape_/shape_); }

    /// Output basic information about given Gaussian RBF.
    template <typename V>
    friend std::ostream& operator<<(std::ostream& os, const Gaussian<V>& m);
};

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_GAUSSIAN_FWD_HPP_
