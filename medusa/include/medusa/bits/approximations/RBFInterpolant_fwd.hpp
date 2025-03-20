#ifndef MEDUSA_BITS_APPROXIMATIONS_RBFINTERPOLANT_FWD_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_RBFINTERPOLANT_FWD_HPP_

/**
 * @file
 * RBFInterpolant class definition.
 *
 * @example test/approximations/RBFInterpolant_test.cpp
 */

#include <medusa/Config.hpp>
#include <Eigen/Core>
#include "Monomials.hpp"

namespace mm {

/**
 * Class representing a RBF Interpolant over a set of nodes of the form
 * @f[
 *    u(p) = \sum_{i=1}^n \alpha_i \phi\left(\left\|\frac{p - p_i}{s}\right\|\right) +
 *           \sum_{i=1}^s \beta_i q_i\left(\frac{p - p_c}{s} \right),
 * @f]
 * where @f$\phi@f$ is the RBF and @f$q_i@f$ are monomials.
 *
 * @note This class does not compute the interpolant, that is done by RBFFD::getApproximant.
 * @sa RBFFD
 * @tparam basis_t Radial basis function used. Must satisfy the @ref rbf-concept.
 * @tparam vec_t Vector type used.
 *
 * Usage example:
 * @snippet test/approximations/RBFInterpolant_test.cpp RBFInterpolant usage example
 */
template <typename RBFType, typename vec_t>
class RBFInterpolant {
  public:
    typedef typename vec_t::scalar_t scalar_t;  ///< Scalar type.
    typedef vec_t vector_t;  ///< Vector type.
    typedef RBFType rbf_t;  ///< Radial basis function type.

  private:
    rbf_t rbf_;  ///< RBF used for interpolation.
    Monomials<vec_t> mon_;  ///< Augmenting monomials.
    vector_t point_;  ///< Center point.
    std::vector<vector_t> support_;  ///< Local scaled stencil points.
    scalar_t scale_;   ///< Scale.
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> coefficients_;  ///< Coefficients (for scaled fn.)

  public:
    /**
     * Construct a RBF interpolant with known coefficients. The interpolant is of the form
     * @f[
     *    u(p) = \sum_{i=1}^n \alpha_i \phi\left(\left\|\frac{p - p_i}{s}\right\|\right) +
     *           \sum_{i=1}^s \beta_i q_i\left(\frac{p - p_c}{s} \right),
     * @f]
     * where @f$\phi@f$ is the RBF and @f$q_i@f$ are monomials.
     * @param rbf Radial basis function @f$\phi@f$, must satisfy the @ref rbf-concept.
     * @param mon Monomial basis @f$q_i@f$.
     * @param point Center point @f$p_c@f$ of the RBF interpolation (used for monomial shift).
     * @param support Nonscaled stencil of the RBFFD approximation.
     * @param scale Scale @f$s@f$ used in RBFFD computation, see @ref scale-concept.
     * @param coefficients The @f$n+s@f$ coefficients @f$\alpha_i, \beta_i@f$ of the approximation.
     */
    RBFInterpolant(const rbf_t& rbf, const Monomials<vec_t>& mon, const vector_t& point,
                   const std::vector<vector_t>& support, scalar_t scale,
                   const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>& coefficients);


    /// Evaluate the interpolant at given `point`.
    scalar_t operator()(const vector_t& point) const;

    /// Evaluate an operator applied to the interpolant at given `point`.
    template <typename operator_t>
    scalar_t operator()(const vector_t& point, const operator_t& op) const;

    /// Get the center point.
    const vector_t& point() const { return point_; }
    /// Get the scale.
    scalar_t scale() const { return scale_; }
    /// Get the coefficient vector.
    const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>& coefficients() const { return coefficients_; }
};


}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_RBFINTERPOLANT_FWD_HPP_
