#ifndef MEDUSA_BITS_APPROXIMATIONS_MONOMIALS_FWD_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_MONOMIALS_FWD_HPP_

/**
 * @file
 * Declaration of Monomial basis.
 *
 * @example test/approximations/Monomials_test.cpp
 */

#include <medusa/Config.hpp>
#include <cmath>
#include <array>
#include <Eigen/Core>
#include <iosfwd>

namespace mm {

template <int dimension>
struct Lap;
template <int dimension>
struct Der1;
template <int dimension>
struct Der2;
template <int dimension>
struct Derivative;

/**
 * A class representing Monomial basis.
 *
 * This class satisfies the @ref basis-concept.
 *
 * Usage example:
 * @snippet approximations/Monomials_test.cpp Monomials usage example
 * @ingroup approximations
 */
template <class vec_t>
class Monomials {
  public:
    typedef vec_t vector_t;   ///< Vector type.
    typedef typename vector_t::scalar_t scalar_t;   ///< Scalar type.
    /// Store dimension of the domain.
    enum { /** Dimensionality of the function domain. */ dim = vec_t::dim };

  private:
    /// A vector describing monomials with the powers of every coordinate
    Eigen::Matrix<int, dim, Eigen::Dynamic> powers_;

    /// Constructs basis from a vector of powers.
    void setFromPowers(const std::vector<std::vector<int>>& powers);

    /// Generate powers for `dim`-d monomials up to `max_order`.
    static std::vector<std::vector<int>> generatePowers(int max_order, int dim);

  public:
    /// Construct empty monomial basis with `size = 0`.
    Monomials() = default;

    /**
     * Construct a basis of all monomials of combined order lower or equal to `order`.
     * @param order Maximal combined order of monomials to be used.
     *  If `order` is `-1` then empty basis is constructed.
     *
     * Example: If you call this with an `order = 2` parameter and Vec2d template parameter,
     * your basis will consist of @f$ \{1, x, y, x^2, xy, y^2\} @f$ (not necessarily in that order).
     */
    Monomials(int order);

    /**
     * @brief Construct monomial basis from monomials with specific powers.
     *
     * @param powers List of lists of size `dim`, each representing powers of a monomial.
     * For example `{{1, 2}, {0, 3}, {2, 0}}` in 2D represents monomials @f$\{x y^2, y^3, x^2\}@f$.
     */
    Monomials(const std::vector<std::vector<int>>& powers) { setFromPowers(powers); }

    /**
     * Construct a tensor basis of monomials ranging from 0 up to `order` (inclusive) in each
     * dimension. If `order` is `-1` then empty basis is constructed.
     *
     * Example: `tensorBasis(2)` in 2D constructs the set
     * @f$\{1, x, x^2, y, yx, yx^2, y^2, y^2x, y^2x^2\}@f$.
     */
    static Monomials<vec_t> tensorBasis(int order);

    /// Return number of monomials in this basis.
    int size() const { return powers_.cols(); }

    /// Get saved monomial powers.
    const Eigen::Matrix<int, dim, Eigen::Dynamic>& powers() const { return powers_; }

    /**
     * Evaluates `index`-th monomial' at `point`.
     *
     * @param index A number in `[0, size())` specifying the index of a monomial to evaluate.
     * @param point Point in which to evaluate the monomial.
     *
     * @return Value of requested monomial at given point.
     * @throws Assertion fails if `index` is out of range or if an invalid derivative is requested.
     */
    scalar_t eval(int index, const vector_t& point,
                  const std::vector<vector_t>& /* support */ = {}) const;

    /**
     * Apply an operator at a given point. For @f$\frac{d}{dx}@f$, the function computes
     * @f$ \frac{d}{dx} p_i(x_s), x_s = \frac{x - c}{s} @f$.
     * @param index A number in `[0, size())` specifying the index of a monomial to evaluate.
     * @param point Translated and scaled point @f$x_s@f$ at which the basis function is evaluated.
     * @param op The differential operator.
     * @param support Translated and scaled values of nodes in the support domain. Unused.
     * @param scale The scaling factor @f$s@f$.
     * @return Value of the operator applied to `index`-th basis function at point `point`.
     * @throws Assertion fails if `index` is out of range.
     * @sa Operator, Lap, Der1, Der2
     */
    template <typename operator_t>
    scalar_t evalOp(int index, const vector_t& point, operator_t op,
                    const std::vector<vector_t>& support = {}, scalar_t scale = 1) const;
    /// Evaluate Laplacian. @sa evalOp
    scalar_t evalOp(int index, const vector_t& point, Lap<dim> op,
                    const std::vector<vector_t>& support = {}, scalar_t scale = 1) const;
    /// Evaluate 1st derivative. @sa evalOp
    scalar_t evalOp(int index, const vector_t& point, Der1<dim> op,
                    const std::vector<vector_t>& support = {}, scalar_t scale = 1) const;
    /// Evaluate 2nd derivative. @sa evalOp
    scalar_t evalOp(int index, const vector_t& point, Der2<dim> op,
                    const std::vector<vector_t>& support = {}, scalar_t scale = 1) const;
    /// Evaluate general derivative. @sa evalOp
    scalar_t evalOp(int index, const vector_t& point, Derivative<dim> op,
                    const std::vector<vector_t>& support = {}, scalar_t scale = 1) const;

    /// Evaluate `index`-th monomial at zero.
    scalar_t evalAt0(int index, const std::vector<vector_t>& /* support */ = {}) const;
    /// Evaluate an operator on `index`-th monomial at zero. @sa eval
    template <typename operator_t>
    scalar_t evalOpAt0(int index, const operator_t& op, const std::vector<vector_t>& support = {},
                       scalar_t scale = 1) const {
        return op.applyAt0(*this, index, support, scale);
    }
    /// Evaluate Laplacian `lap` at zero. @sa evalOpAt0
    scalar_t evalOpAt0(int index, const Lap<dim>& lap, const std::vector<vector_t>& = {},
                       scalar_t scale = 1) const;
    /// Evaluate first derivative operator `der1` at zero. @sa evalOpAt0
    scalar_t evalOpAt0(int index, const Der1<dim>& der1, const std::vector<vector_t>& = {},
                       scalar_t scale = 1) const;
    /// Evaluate second derivative operator `der2` at zero. @sa evalOpAt0
    scalar_t evalOpAt0(int index, const Der2<dim>& der2, const std::vector<vector_t>& = {},
                       scalar_t scale = 1) const;

    /// Evaluate general derivative operator at zero. @sa evalOpAt0
    scalar_t evalOpAt0(int index, const Derivative<dim>& der, const std::vector<vector_t>& = {},
                       scalar_t scale = 1) const;

    /// Output basic info about given Monomial basis.
    template <class V>
    friend std::ostream& operator<<(std::ostream& os, const Monomials<V>& m);
};

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_MONOMIALS_FWD_HPP_
