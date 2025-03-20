#ifndef MEDUSA_BITS_APPROXIMATIONS_OPERATORS_FWD_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_OPERATORS_FWD_HPP_

/**
 * @file
 * File defining basic differential operator families.
 *
 * @example test/approximations/operators_test.cpp
 */

#include <medusa/Config.hpp>
#include <array>
#include <vector>
#include <Eigen/Core>
#include "Monomials_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>

namespace mm {

/**
 * Base class for a differential operator. Also doubles as a one-element family.
 * @tparam Derived Derived class inheriting, used for CRTP.
 * @ingroup approximations
 *
 * Example of custom operator usage. Definition:
 * @snippet approximations/operators_test.cpp Biharmonic def
 * Usage:
 * @snippet approximations/operators_test.cpp Biharmonic usage
 *
 * See also custom_operators example in the examples suite.
 */
template <typename Derived>
struct Operator {
    typedef Derived operator_t;  ///< Type of the contained operators.
    static constexpr int size() { return 1; }   ///< Size of a single-operator family is always 1.
    static constexpr int index(Derived) { return 0; }   ///< Index is always 0.
    std::string name() const { return Derived::name(); }   ///< CRTP for human-readable name.
    static std::string type_name() { return Derived::type_name(); }  ///< Human-readable type name.

    /**
     * Return an iterable of iterators represented by a family.
     * As this is a single operator, the family consists of itself only.
     */
    std::array<Derived, 1> operators() const { return {*static_cast<const Derived*>(this)}; }

    /**
     * Apply this operator to a given basis function. This function should be specialized
     * (overloaded) for more concrete values of `basis_t` in concrete classes inheriting form
     * this one. It is intended that this function is never called but only serves as a default
     * for demonstration purposes and better error messages.
     * @tparam basis_t Basis type.
     * @param basis The basis @f$b@f$ to apply the operator to.
     * @param index Which element of the basis to apply the operator to.
     * @param point At which aldeary shifted and translated point to apply the operator.
     * @param support Local support domain (already shifted and scaled). These are the values
     * @f$(x-c)/s@f$ in the equation below.
     * @param scale Local scaling factor @f$s@f$
     * @return Let this operator be called @f$L@f$, the `index`-th basis function @f$b_i@f$.
     * The function returns @f$L b_i((x-c)/s)@f$. Not that the operator @f$L@f$ is usually
     * something like "derivative wrt. x" and factors of @f$s@f$ appear in the result.
     * @sa RBFBasis::evalOp, Monomials::evalOp, applyAt0
     */
    template <typename basis_t>
    typename basis_t::scalar_t apply(
            const basis_t& basis, int index, typename basis_t::vector_t point,
            const std::vector<typename basis_t::vector_t>& support,
            typename basis_t::scalar_t scale) const;

    /**
     * Like @ref apply, but with point equal to 0.
     * @sa apply, RBFBasis::evalOpAt0, Monomials::evalOpAt0
     */
    template <typename basis_t>
    typename basis_t::scalar_t applyAt0(
            const basis_t& basis, int index,
            const std::vector<typename basis_t::vector_t>& support,
            typename basis_t::scalar_t scale) const;
};
/// Output info about given operator.
template <typename Derived>
std::ostream& operator<<(std::ostream& os, const Operator<Derived>& op) { return os << op.name(); }

/// Represents the Laplacian operator. @ingroup approximations
template <int dimension>
struct Lap : public Operator<Lap<dimension>> {
    /// Store dimension of the domain.
    enum { /** Dimensionality of the domain. */ dim = dimension };
    static std::string name() { return format("Laplacian<%d>", dim); }  ///< Human-readable name.
    static std::string type_name() { return name(); }  ///< Human readable type name.
};

/**
 * Represents a family of all first derivatives.
 * @ingroup approximations
 */
template <int dimension>
struct Der1s {
    /// Store dimension of the domain.
    enum { /** Dimensionality of the domain. */ dim = dimension };
    typedef Der1<dim> operator_t;  ///< Type of this family's members.
    static int size() { return dim; }  ///< Size of this family is `dim`.
    static std::string name() { return format("Der1s<%d>", dim); }  ///< Human-readable name.
    static std::string type_name() { return name(); }  ///< Human readable type name.

    /// Return an iterable of iterators represented by a family.
    static std::array<Der1<dimension>, dimension> operators();

    /// Get index of operator `d` in the family.
    static int index(Der1<dim> d) { return d.var; }
};
/// Output info about Der1s<d>
template <int d>
std::ostream& operator<<(std::ostream& os, const Der1s<d>& op) { return os << op.name(); }

/**
 * Represents a family of all second derivatives.
 * @ingroup approximations
 */
template <int dimension>
struct Der2s {
    /// Store dimension of the domain.
    enum { /** Dimensionality of the domain. */ dim = dimension };
    typedef Der2<dim> operator_t;  ///< Type of this family's members.
    static int size() { return dim * (dim + 1) / 2; }  ///< Size of this family is `dim*(dim-1)/2`.
    static std::string name() { return format("Der2s<%d>", dim); }  ///< Human-readable name.
    static std::string type_name() { return name(); }  ///< Human readable type name.

    /// Return an iterable of iterators represented by a family.
    static std::array<Der2<dimension>, dimension * (dimension + 1) / 2> operators();

    /// Get index of operator `d` in the family.
    static int index(Der2<dim> d) { return (d.var2 == 2 ? d.var1 + 1 : d.var1) + d.var2; }
};
/// Output info about Der2s<d>
template <int d>
std::ostream& operator<<(std::ostream& os, const Der2s<d>& op) { return os << op.name(); }


/// Represents a first derivative wrt. `var`. @ingroup approximations
template <int dimension>
struct Der1 : public Operator<Der1<dimension>> {
    /// Store dimension of the domain.
    enum { /** Dimensionality of the domain. */ dim = dimension };
    int var;  ///< Index representing derived variable (x = 0, y = 1, z = 2)

    /// Default constructor for array preallocation purposes.
    Der1() : var(-1) {}

    /**
     * Initialization of specific 1st derivative operator.
     * @param var Index representing derived variable (x = 0, y = 1, z = 2).
     */
    Der1(int var);

    /// Human readable name.
    std::string name() const;
    /// Human readable type name.
    static std::string type_name() { return format("Der1<%d>", dim); }
};


/// Represents a second derivative wrt. `var1` and `var2`. @ingroup approximations
template <int dimension>
struct Der2 : public Operator<Der2<dimension>> {
    /// Store dimension of the domain.
    enum { /** Dimensionality of the domain. */ dim = dimension };
    int var1,  ///< Lower index representing derived variable (x = 0, y = 1, z = 2, ...)
        var2;  ///< Higher index representing derived variable (x = 0, y = 1, z = 2, ...)

    /// Default constructor for array preallocation purposes.
    Der2() : var1(-1), var2(-1) {}

    /**
     * Initialization of specific 2nd derivative operator.
     * @param var1 First variable with respect to which the derivative is taken.
     * variable (x = 0, y = 1, z = 2, ...).
     * @param var2 Second variable with respect to which the derivative is taken.
     * (Should have higher index and the first variable.)
     * variable (x = 0, y = 1, z = 2, ...).
     */
    Der2(int var1, int var2);

    /**
     * Initialization of specific 2nd derivative operator.
     * @param var Index representing two times derived variable (x = 0, y = 1, z = 2, ...).
     */
    Der2(int var) : Der2(var, var) {}

    /// Human readable name.
    std::string name() const;
    /// Human readable type name.
    static std::string type_name() { return format("Der2<%d>", dim); }
};


/**
 * Represents a general derivative @f$\frac{\partial^{|\alpha|}}{\partial x^\alpha}@f$.
 * @ingroup approximations
 */
template <int dimension>
struct Derivative : public Operator<Derivative<dimension>> {
    /// Store dimension of the domain.
    enum { /** Dimensionality of the domain. */ dim = dimension };
    std::array<int, dim> orders;  ///< Values @f$\alpha@f$.

    /// Parameter `orders` represent the derivative orders @f$\alpha@f$.
    Derivative(const std::array<int, dim>& orders) : orders(orders) {}

    /// Human readable name.
    std::string name() { return format("Derivative<%d> wrt. %s", dim, orders); }
    /// Human readable type name.
    static std::string type_name() { return format("Derivative<%d>", dim); }
};

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_OPERATORS_FWD_HPP_
