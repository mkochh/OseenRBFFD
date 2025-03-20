#ifndef MEDUSA_BITS_APPROXIMATIONS_WLSAPPROXIMANT_FWD_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_WLSAPPROXIMANT_FWD_HPP_

/**
 * @file
 * Declaration of the function, obtained by WLS approximation.
 *
 * @example test/approximations/WLSApproximant_test.cpp
 */

#include <medusa/Config.hpp>
#include <Eigen/Core>

namespace mm {

/**
 * Class representing the function that is a WLS approximant using some
 * basis function over some points. These must satisfy the @ref basis-concept.
 * @note This class does not compute the approximant, that is done by WLS::getApproximant.
 * @sa WLS
 * @tparam basis_t Basis function type.
 *
 * Usage example:
 * @snippet test/approximations/WLSApproximant_test.cpp WLSApproximant usage example
 */
template <typename basis_t>
class WLSApproximant {
  public:
    typedef typename basis_t::scalar_t scalar_t;  ///< Scalar type.
    typedef typename basis_t::vector_t vector_t;  ///< Vector type.

  private:
    basis_t basis_;  ///< Basis.
    vector_t point_;  ///< Center point.
    std::vector<vector_t> support_;  ///< Local scaled stencil points.
    scalar_t scale_;   ///< Scale.
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> coefficients_;  ///< Coefficients (for scaled fn.)
    scalar_t residual_;   ///< Store residual of the approximation.

  public:
    /**
     * Construct a WLS approximant with known coefficients. The approximant is of the form
     * @f[ u(p) = \sum_{i=1}^m \alpha_i b_i\left(\frac{p - p_c}{s}\right), @f]
     * where @f$b_i@f$ may depend on the stencil nodes (e.g. RBF-s).
     * @param basis Basis functions @f$b_i@f$ to be used, must satisfy the @ref basis-concept.
     * @param point Center point @f$p_c@f$ of the WLS approximation.
     * @param support Nonscaled stencil of the WLS approximation.
     * @param scale Scale @f$s@f$ used in WLS computation, see @ref scale-concept.
     * @param coefficients Coefficients @f$\alpha_i@f$ of the approximation.
     * @param residual The resudial of the weighted least squares approximation. If 0, this
     * function is an interpolant.
     */
    WLSApproximant(const basis_t& basis, const vector_t& point,
                   const std::vector<vector_t>& support, scalar_t scale,
                   const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>& coefficients,
                   scalar_t residual = NaN);

    /// Evaluate the approximant at given `point`.
    scalar_t operator()(const vector_t& point) const;

    /// Evaluate an operator applied to approximant at given `point`.
    template <typename operator_t>
    scalar_t operator()(const vector_t& point, const operator_t& op) const;

    /// Get the basis functions.
    const basis_t& basis() const { return basis_; }
    /// Get the center point.
    const vector_t& point() const { return point_; }
    /// Get the scale.
    scalar_t scale() const { return scale_; }
    /// Get the coefficient vector.
    const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>& coefficients() const { return coefficients_; }
    /// Get the residual.
    scalar_t residual() const { return residual_; }
};

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_WLSAPPROXIMANT_FWD_HPP_
