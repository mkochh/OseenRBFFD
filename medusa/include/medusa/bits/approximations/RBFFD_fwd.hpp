#ifndef MEDUSA_BITS_APPROXIMATIONS_RBFFD_FWD_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_RBFFD_FWD_HPP_

/**
 * @file
 * Declaration of Radial Basis Function Finite Difference approximation.
 *
 * @example test/approximations/RBFFD_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/types/Range_fwd.hpp>
#include <Eigen/Core>
#include "RBFBasis_fwd.hpp"
#include "Monomials_fwd.hpp"
#include "ScaleFunction.hpp"
#include "RBFInterpolant_fwd.hpp"

namespace mm {

/**
 * Computes a RBF-FD approximation of given operator over local neighbourhood.
 *
 * @tparam RBFType Type of RBF used in approximation. Must satisfy the @ref rbf-concept.
 * @tparam vec_t Vector type used in calculations.
 * @tparam scale_t Scale function to be used in approximation. Must satisfy the @ref scale-concept.
 * @tparam solver_t Which linear solver to use in weight computations. Must satisfy the
 * @ref linsolve-concept. Default: LU with partial pivoting.
 *
 * This class satisfies the @ref approx-concept.
 *
 * Usage example:
 * @snippet approximations/RBFFD_test.cpp RBFFD usage example
 * @ingroup approximations
 */
template <class RBFType, class vec_t, class scale_t = NoScale, class solver_t =
    Eigen::PartialPivLU<Eigen::Matrix<typename vec_t::scalar_t, Eigen::Dynamic, Eigen::Dynamic>>>
class RBFFD {
    static_assert(std::is_same<typename vec_t::scalar_t, typename RBFType::scalar_t>::value,
                  "Basis and underlying RBF must have the same scalar type.");
  public:
    typedef RBFType rbf_t;   ///< Radial basis function type.
    typedef vec_t vector_t;   ///< Vector type.
    typedef typename vector_t::scalar_t scalar_t;   ///< Scalar type.
    /// Store dimension of the domain.
    enum { /** Dimensionality of the function domain. */ dim = vec_t::dim };

  private:
    rbf_t rbf_;  ///< RBF used in the approximation.
    Monomials<vec_t> mon_;  ///< Monomials used in the approximation.
    solver_t solver_;  ///< Linear solver used in the approximation.

    vector_t point_;  ///< Saved point.
    scalar_t scale_;  ///< Saved scale.
    Range<vec_t> support_;  ///< Saved support scaled to local coordinate system.

  public:
    /// Construct a RBFFD engine from given RBF using no monomial augmentation.
    explicit RBFFD(rbf_t rbf) : rbf_(rbf), mon_(), point_(NaN) {}

    /// Construct a RBFFD engine from given RBF with additional monomial augmentation.
    RBFFD(rbf_t rbf, const Monomials<vec_t>& mon) : rbf_(rbf), mon_(mon), point_(NaN) {}

    /// Return the RBF used in this approximation.
    const rbf_t& rbf() const { return rbf_; }

    /// Return the augmenting monomial basis.
    const Monomials<vec_t>& monomials() const { return mon_; }

    /// Returns minimal support size required for this approximation to be performed.
    int requiredSupportSize() const { return mon_.size(); }

    /**
     * Setup this engine to approximate operators at given point over given local neighbourhood.
     * @param point At which point to construct the approximation.
     * @param support Coordinates of local support points.
     */
    void compute(const vector_t& point, const std::vector<vector_t>& support);

    /// Return shape function (stencil weights) for value approximation.
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> getShape() const;

    /// Return shape function (stencil weights) for given derivative operator.
    template <class operator_t>
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> getShape(const operator_t& op) const;

    /**
     * Returns weights for approximation with a given right hand side, which is expected to
     * contains the values of an operator applied to your basis functions (both RBFs and monomials).
     * This method is useful if values of the operator are computed differently to the procedure
     * used by @ref getShape.
     * @param Lb The right hand side of the approximation system representing an values of
     * an operator applied to basis functions (both RBFs and monomials) at the point where
     * @ref compute was called.
     * @sa getShape
     */
    template <typename Derived>
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> getShapeFromRhs(
            const Eigen::MatrixBase<Derived>& Lb) const;

  private:
    /// Returns the augmented collocation matrix.
    Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> getMatrix(
            const std::vector<vector_t>& local_coordinates) const;

  public:
    /**
     * Return the center point of the approximation.
     * Only initialized after the first call to @ref compute.
     */
    vec_t center() const { return point_; }

    /// Return the scaling factor. Only initialized after the first call to @ref compute.
    scalar_t scale() const { return scale_; }

    /**
     * Returns the local (translated and scaled) coordinates of the support nodes.
     * Only initialized after the first call to @ref compute.
     */
    Range<vec_t> localCoordinates() const { return support_; }

    /**
     * Returns the augmented collocation matrix.
     * Only initialized after the first call to @ref compute.
     */
    Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> getMatrix() const {
        return getMatrix(support_);
    }

    /**
     * Returns the stored decomposition of the augmented collocation matrix.
     * Only initialized after the first call to @ref compute.
     */
    const solver_t& solver() const { return solver_; }

    /**
     * Construct a RBF interpolant of form
     * @f[
     *    u(p) = \sum_{i=1}^n \alpha_i \phi\left(\left\|\frac{p - p_i}{s}\right\|\right) +
     *           \sum_{i=1}^s \beta_i q_i\left(\frac{p - p_c}{s} \right),
     * @f]
     * where @f$\phi@f$ is the RBF and @f$q_i@f$ are monomials. Additionally, @f$s@f$
     * constraints are imposed for @f$\beta_i@f$, @f$\sum_{i=1}^n q_j(p_i) = 0@f$, for
     * @f$j = 1, \ldots, s@f$.
     * @param point Center point @f$p@f$.
     * @param support Support values @f$p_i@f$.
     * @param values Function values @f$u_i@f$ at @f$p_i@f$.
     * @return An object representing the interpolant.
     */
    template <typename Derived>
    RBFInterpolant<rbf_t, vec_t> getApproximant(const vector_t& point,
            const std::vector<vector_t>& support, const Eigen::MatrixBase<Derived>& values) const;

    /// Output basic info about given approximation engine.
    template <class R, class V, class S, class L>
    friend std::ostream& operator<<(std::ostream& os, const RBFFD<R, V, S, L>& e);
};

/// Output basic info about given approximation engine.
template <class R, class V, class S, class L>
std::ostream& operator<<(std::ostream& os, const RBFFD<R, V, S, L>& e) {
    os << "RBFFD approximation engine:\n"
       << "    dimension: " << e.dim << '\n';
    if (e.point_[0] != e.point_[0]) {
        os << "    last point: not used yet\n";
    } else {
        os << "    last point: " << e.point_ << '\n';
    }
    os << "    RBF: " << e.rbf_ << '\n'
       << "    Augmentation: " << e.mon_ << '\n'
       << "    scale: " << S();
    return os;
}

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_RBFFD_FWD_HPP_
