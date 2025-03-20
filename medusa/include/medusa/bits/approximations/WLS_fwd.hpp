#ifndef MEDUSA_BITS_APPROXIMATIONS_WLS_FWD_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_WLS_FWD_HPP_

/**
 * @file
 * Declaration of weighted least squares approximation.
 *
 * @example test/approximations/WLS_test.cpp WLS
 */

#include <medusa/Config.hpp>
#include <medusa/bits/types/Range_fwd.hpp>
#include <Eigen/Core>
#include <type_traits>
#include <iosfwd>
#include "WLSApproximant_fwd.hpp"
#include "Operators_fwd.hpp"

namespace mm {

template <typename scalar_t, int QRPreconditioner>
class JacobiSVDWrapper;

template <typename vec_t>
class NoWeight;

class NoScale;

/**
 * A class for generating approximations using Weighted Least Squares over local neighborhoods.
 *
 * @tparam basis_t Type of basis used in approximation. Must satisfy the @ref basis-concept.
 * @tparam weight_t Type of weight function used in approximation. Must satisfy the
 * @ref weight-concept.
 * @tparam scale_t Scale function to be used in approximation. Must satisfy the @ref scale-concept.
 * @tparam solver_t Which linear solver to use in weight computations. Must satisfy the
 * @ref linsolve-concept. Default: Jacobi SVD.
 *
 * This class satisfies the @ref approx-concept.
 *
 * Usage example:
 * @snippet approximations/WLS_test.cpp WLS usage example
 * @ingroup approximations
 */
template <class basis_t, class weight_t = NoWeight<typename basis_t::vector_t>,
          class scale_t = NoScale, class solver_t =
               JacobiSVDWrapper<typename basis_t::scalar_t,
                                Eigen::ColPivHouseholderQRPreconditioner>>
class WLS {
  public:
    static_assert(std::is_same<typename basis_t::vector_t, typename weight_t::vector_t>::value,
                  "Basis and weight functions must have the same vector type.");

    /// Numeric scalar data type, e.g.\ double.
    typedef typename basis_t::scalar_t scalar_t;
    /// Vector data type.
    typedef typename basis_t::vector_t vector_t;
    /// Store dimension of the domain.
    enum { /** Dimensionality of the domain. */ dim = vector_t::dim };
    /// Eigen matrix data type.
    typedef typename Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> ei_matrix_t;

  private:
    basis_t basis_;  ///< Basis functions.
    weight_t weight_;  ///< Weight function.

    Range<vector_t> local_coordinates_;  ///< Support points after translation and scaling.
    vector_t point_;  ///< Saved center point.
    scalar_t scale_;  ///< Saved scale.
    solver_t solver_;  ///< Saved decomposition.

  public:
    /// Construct WLS with given basis and default constructed weight.
    explicit WLS(const basis_t& basis) : basis_(basis), scale_(NaN) {}
    /**
     * Hold information about the type of approximation we want to have for the domain.
     * @param basis Basis to use in the algorithm.
     * @param weight Weight function to use.
     */
    WLS(const basis_t& basis, const weight_t& weight) :
            basis_(basis), weight_(weight), scale_(NaN) {}

  private:
    /// Returns the collocation matrix.
    ei_matrix_t getMatrix(const std::vector<vector_t>& local_coordinates) const;

  public:
    /**
     * Computes, decomposes and stores matrix for shape calculation at given point
     * with given neighbours. Given basis and weight function are used
     * and evaluated in the local coordinate system with origin in `point` and
     * scales by @ref scale_. Decomposition given by @ref solver_ is used.
     * This function must be called before any calls to @ref getShape are made for this
     * point configuration.
     * @param point Point around which to prepare for shape computation.
     * @param support Neighbouring points where function values are expected.
     */
    void compute(const vector_t& point, const std::vector<vector_t>& support);

    /**
     * Returns weights for approximation of function value at point `point`, knowing
     * values in its `support`, given by a previous call to compute.
     * @return Vector of nodal weights @f$\chi@f$, such that
     * @f$u(p) \approx \sum_{i=1}^n \chi_i u(p_i)@f$ where @f$p@f$ represents `point` and
     * @f$p_i@f$ are its support nodes.
     */
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> getShape() const;

    /**
     * Returns weights for approximation of differentiation operator, given by `der`
     * at point `point`, knowing values in its `support`, given by a previous call to compute.
     * @return Vector of nodal weights @f$\chi@f$, such that
     * @f$(\mathcal{L}u)(p) \approx \sum_{i=1}^n \chi_i u(p_i)@f$ where @f$p@f$ represents
     * `point`, @f$p_i@f$ are its support nodes and @f$\mathcal{L}@f$ is given by `der` as
     * @f$\mathcal{L} = \dpar{^{\sum_{i=1}^d m_i}}{x_1^{m_1} \cdots \partial x_d^{m_d}}@f$,
     * where @f$m_i@f$ are elements of `der`.
     * @param op Operator to be approximated.
     */
    template <class operator_t>
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> getShape(const operator_t& op) const;

    /**
     * Returns weights for approximation with a given right hand side, which is expected to
     * contains the values of an operator applied to your basis functions. This method
     * is useful if values of the operator are computed differently to the procedure used by
     * @ref getShape.
     * @param Lb The right hand side of the approximation system representing an values of
     * an operator applied to basis functions at the point where @ref compute was called.
     * @sa getShape
     */
    template <typename Derived>
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> getShapeFromRhs(
            const Eigen::MatrixBase<Derived>& Lb) const;

    /// Returns minimal support size required for this approximation to be performed.
    int requiredSupportSize() const { return basis_.size(); }
    /// Returns basis used in the computation.
    const basis_t& basis() const { return basis_; }
    /// Returns weight function used in the computation.
    const weight_t& weight() const { return weight_; }

    /**
     * Return the center point of the approximation.
     * Only initialized after the first call to @ref compute.
     */
    vector_t center() const { return point_; }

    /// Return the scaling factor. Only initialized after the first call to @ref compute.
    scalar_t scale() const { return scale_; }

    /**
     * Returns the local (translated and scaled) coordinates of the support nodes.
     * Only initialized after the first call to @ref compute.
     */
    Range<vector_t> localCoordinates() const { return local_coordinates_; }

    /**
     * Returns the matrix used in the approximation.
     * Only initialized after the first call to @ref compute.
     */
    ei_matrix_t getMatrix() const { return getMatrix(local_coordinates_); }

    /// Returns the current support weights. Only initialized after the first call to @ref compute.
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> getWeights() const;

    /**
     * Read access to the linear solver used for matrix decomposition.
     * Only initialized after the first call to @ref compute.
     */
    const solver_t& solver() const { return solver_; }

    /**
     * Read-write access to the linear solver used for matrix decomposition.
     * Can be used to set additional solver parameters.
     * Only initialized after the first call to @ref compute.
     */
    solver_t& solver() { return solver_; }

    /**
     * Construct a WLS approximant anchored at `point` with given values in given support nodes.
     * The approximation function @f$\hat{u}(x) =
     * \sum_{j=1}^m \alpha_j b_j\left(\frac{x - p}{s}\right)@f$ minimises
     * @f[ J = \sum_{i=1}^n w\left( \frac{p_i-p}{s}\right)^2 (u_i - \hat{u}(p_i))^2 @f]
     * @param point Center point @f$p@f$.
     * @param support Support values @f$p_i@f$.
     * @param values Function values @f$u_i@f$ at @f$p_i@f$.
     * @return An object representing the approximant.
     */
    template <typename Derived>
    WLSApproximant<basis_t> getApproximant(const vector_t& point,
            const std::vector<vector_t>& support, const Eigen::MatrixBase<Derived>& values) const;

    /// Output basic info about given approximation engine.
    template <class B, class W, class S, class L>
    friend std::ostream& operator<<(std::ostream& os, const WLS<B, W, S, L>& wls);
};

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_WLS_FWD_HPP_
