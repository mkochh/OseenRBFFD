#ifndef MEDUSA_BITS_APPROXIMATIONS_RBFBASIS_FWD_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_RBFBASIS_FWD_HPP_

/**
 * @file
 * Declaration of RBF basis.
 *
 * @example test/approximations/RBFBasis_test.cpp
 */

#include <medusa/Config.hpp>
#include <iosfwd>
#include <array>
#include "Operators_fwd.hpp"

namespace mm {

/**
 * Represents a basis of Radial Basis Functions over a local neighbourhood.
 * @tparam RBFType Type of the RBF used. Must satisfy the @ref rbf-concept.
 * @tparam vec_t Vector type used in calculations.
 *
 * This class satisfies the @ref basis-concept.
 *
 * Usage example:
 * @snippet approximations/RBFBasis_test.cpp RBF basis usage example
 * @ingroup approximations
 */
template <class RBFType, class vec_t>
class RBFBasis {
    static_assert(std::is_same<typename vec_t::scalar_t, typename RBFType::scalar_t>::value,
                  "Basis and underlying RBF must have the same scalar type.");
  public:
    typedef RBFType rbf_t;   ///< Radial basis function type.
    typedef vec_t vector_t;   ///< Vector type.
    typedef typename vector_t::scalar_t scalar_t;   ///< Scalar type.
    /// Store dimension of the domain.
    enum { /** Dimensionality of the function domain. */ dim = vec_t::dim };

  private:
    int size_;  ///< Basis size.
    rbf_t rbf_;  ///< RBF function.

  public:
    /// Construct a basis of `basis_size` with default construction of RBF.
    explicit RBFBasis(int basis_size) : size_(basis_size), rbf_() {}
    /// Construct a basis of `basis_size` with given RBF.
    RBFBasis(int basis_size, const rbf_t& rbf) : size_(basis_size), rbf_(rbf) {}

    /// Returns basis size.
    int size() const { return size_; }

    /// Returns underlying RBF object.
    const rbf_t& rbf() const { return rbf_; }

    /// Returns modifiable underlying RBF object.
    rbf_t& rbf() { return rbf_; }


    /**
     * Evaluates `index`-th RBF's at `point`.
     *
     * @param index A number in `[0, size())` specifying the index of RBF to evaluate.
     * @param point Point in which to evaluate the RBF.
     * @param support Points in local neighbourhood used as centers for RBFs.
     *
     * @return Value of requested monomial at given point.
     * @throws Assertion fails if `index` is out of range or if an invalid derivative is requested.
     *
     * @note Derivatives only up to combined order of 2 are currently supported.
     */
    scalar_t eval(int index, const vector_t& point, const std::vector<vector_t>& support) const;

    /**
     * Apply an operator at a given point. For @f$\frac{d}{dx}@f$, the function computes
     * @f$ \frac{d}{dx} p_i(x_s), x_s = \frac{x - c}{s} @f$.
     * @param index A number in `[0, size())` specifying the index of a monomial to evaluate.
     * @param point Translated and scaled point @f$x_s@f$ at which the basis function is evaluated.
     * @param op The differential operator.
     * @param support Translated and scaled values of nodes in the support domain.
     * @param scale The scaling factor @f$s@f$.
     * @return Value of the operator applied to `index`-th basis function at point `point`.
     * @throws Assertion fails if `index` is out of range.
     * @sa Operator, Lap, Der1, Der2
     */
    template <typename operator_t>
    scalar_t evalOp(int index, const vector_t& point, operator_t op,
                    const std::vector<vector_t>& support, scalar_t scale = 1) const;
    /// Evaluate Laplacian. @sa evalOp
    scalar_t evalOp(int index, const vector_t& point, Lap<dim> op,
                    const std::vector<vector_t>& support, scalar_t scale = 1) const;
    /// Evaluate 1st derivative. @sa evalOp
    scalar_t evalOp(int index, const vector_t& point, Der1<dim> op,
                    const std::vector<vector_t>& support, scalar_t scale = 1) const;
    /// Evaluate 2nd derivative. @sa evalOp
    scalar_t evalOp(int index, const vector_t& point, Der2<dim> op,
                    const std::vector<vector_t>& support, scalar_t scale = 1) const;

    /// Evaluate `index`-th function at zero.
    scalar_t evalAt0(int index, const std::vector<vector_t>& support) const;

    /**
     * Evaluate general operator `op` at zero.
     *
     * @param index A number in `[0, size())` specifying the index of RBF to evaluate.
     * @param op Evaluated operator.

     * @param support Points in local neighbourhood used as centers for RBFs.
     * @param scale Scaling factor.
     *
     * @return Value of requested monomial at given point.
     * @throws Assertion fails if `index` is out of range or if an invalid derivative is requested.
     */
    template <typename operator_t>
    scalar_t evalOpAt0(int index, const operator_t& op, const std::vector<vector_t>& support,
                       scalar_t scale = 1) const {
        return op.applyAt0(*this, index, support, scale);
    }
    /// Evaluate Laplacian `lap` at zero.
    scalar_t evalOpAt0(int index, const Lap<dim>& lap, const std::vector<vector_t>& support,
                       scalar_t scale = 1) const;
    /// Evaluate first derivative operator `der1` at zero.
    scalar_t evalOpAt0(int index, const Der1<dim>& der1, const std::vector<vector_t>& support,
                       scalar_t scale = 1) const;
    /// Evaluate second derivative operator `der2` at zero.
    scalar_t evalOpAt0(int index, const Der2<dim>& der2, const std::vector<vector_t>& support,
                       scalar_t scale = 1) const;

    /// Output basic info about given basis.
    template <typename V, typename R>
    friend std::ostream& operator<<(std::ostream& os, const RBFBasis<V, R>& m);
};


// Convenience typedefs

template <typename S> class Gaussian;
/// RBF basis using Gaussian RBF. Defined for convenience.
template <typename V> using Gaussians = RBFBasis<Gaussian<typename V::scalar_t>, V>;

template <typename S> class Multiquadric;
/// RBF basis using Multiquadric RBF. Defined for convenience.
template <typename V> using MQs = RBFBasis<Multiquadric<typename V::scalar_t>, V>;

template <typename S> class InverseMultiquadric;
/// RBF basis using InverseMultiquadric RBF. Defined for convenience.
template <typename V> using IMQs = RBFBasis<InverseMultiquadric<typename V::scalar_t>, V>;

template <typename S, int K> class Polyharmonic;
/// RBF basis using Polyharmonic RBF. Defined for convenience.
template <typename V, int K = -1> using PHs = RBFBasis<Polyharmonic<typename V::scalar_t, K>, V>;

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_RBFBASIS_FWD_HPP_
