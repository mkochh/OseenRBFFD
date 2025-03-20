#ifndef MEDUSA_BITS_OPERATORS_EXPLICITOPERATORS_FWD_HPP_
#define MEDUSA_BITS_OPERATORS_EXPLICITOPERATORS_FWD_HPP_

/**
 * @file
 * Explicit operators declarations.
 *
 * @example test/operators/ExplicitOperators_test.cpp
 */

#include <medusa/Config.hpp>
#include <Eigen/Core>
#include "UniformShapeStorage_fwd.hpp"
#include "RaggedShapeStorage_fwd.hpp"
#include <medusa/bits/types/Vec.hpp>

namespace mm {

/**
 * A class for evaluating typical operators needed in spatial discretization.
 * This class represents scalar differential operators, discretely approximated
 * with shape functions that are stored in given shape storage. These shape functions
 * are used to construct approximations of common operators such as gradient, Laplacian
 * and coordinate derivatives that can be applied to scalar fields at given points.
 *
 * @tparam shape_storage_type Storage class for shapes. @ref ss-concept.
 *
 * Usage example:
 * @snippet operators/ExplicitOperators_test.cpp Explicit operators usage example
 * @ingroup operators
 *
 * @sa ExplicitVectorOperators
 */
template <class shape_storage_type>
class ExplicitOperators {
  public:
    typedef shape_storage_type shape_storage_t;   ///< Type of shape storage.
    typedef typename shape_storage_t::vector_t vector_t;  ///< Vector type.
    typedef typename shape_storage_t::scalar_t scalar_t;  ///< Scalar type.
    /// Store dimension of the domain.
    enum { /** Dimensionality of the function domain. */ dim = shape_storage_t::dim };

  private:
    /// Pointer to shape storage, but name is shortened for readability.
    const shape_storage_t* ss;

  public:
    /// Construct empty explicit operators.
    ExplicitOperators() : ss(nullptr) {}
    /**
     * Construct explicit operators over given shape storage.
     * The pointer to given storage is stored as a non-owning pointer. It is the user's
     * responsibility that given storage is valid throughout the use use of this class.
     *
     * This class is usually constructed directly from shape storage using the `explicitOperators()`
     * member function.
     */
    explicit ExplicitOperators(const shape_storage_t& ss) : ss(&ss) {}

    /**
     * Sets a new shape storage from which this operators are generated. As in the constructor,
     * the pointer to given storage is stored as a non-owning pointer. It is the user's
     * responsibility that given storage is valid throughout the use use of this class.
     * @param shape_storage_ New shape storage.
     */
    void setShapes(const shape_storage_t& shape_storage_) { ss = &shape_storage_; }

    /// Returns `true` if operators have a non-null pointer to storage and `false` otherwise.
    bool hasShapes() const { return ss != nullptr; }

    /**
     * Returns an approximation of Laplacian of field `u` in `node`-th node.
     * @returns The Laplacian of `u`, e.g.\ when `dim = 3`:
     * @f$
     * \dpar{^2 u}{ x^2} +
     * \dpar{^2 u}{ y^2} +
     * \dpar{^2 u}{ z^2}
     * @f$, where derivatives are approximated using shape functions from @ref ss.
     */
    template <typename scalar_field_t>
    typename scalar_type<scalar_field_t>::type lap(const scalar_field_t& u, int node) const {
        return apply<Lap<dim>, scalar_field_t>(u, node);
    }

    /**
     * Returns an approximation of gradient of field `u` in `node`-th node.
     * @return Vector containing the components of the gradient, e.g.\ when `dim = 3`:
     * @f$ [
     *   \dpar{ u}{ x},
     *   \dpar{ u}{ y},
     *   \dpar{ u}{ z}
     * ]^T@f$, where derivatives are approximated using shape functions from @ref ss.
     */
    template <typename scalar_field_t>
    Vec<typename scalar_type<scalar_field_t>::type, ExplicitOperators<shape_storage_type>::dim>
    grad(const scalar_field_t& u, int node) const;

    /**
     * Calculates a new field value `a` for `u` at `node`, such that the Neumann boundary
     * condition `du/dn(node) = val` holds.
     * Denote given index `node` by @f$i@f$ and `val` by @f$v@f$.
     * The function returns such @f$a@f$, that by setting @f$ u(i) = a@f$ we satisfy the condition
     * @f$ \dpar{ u}{ \vec n}(i) = v@f$. This is done
     * using the following formula:
     * @f[ a = \frac{v - \sum\limits_{i=1}^{s} \sum\limits_{d=1}^{dim}
     *   n_d \chi_{\dpar{}{ x_d}}(i)\ u(i)}
     * {\sum\limits_{d=1}^{dim} n_d \chi_{\dpar{}{ x_d}}(0)},
     * @f]
     * where @f$\chi_{\dpar{}{ x_d}}@f$ represents the formula for
     * the shape function of first derivative wrt. @f$d@f$.
     *
     * @note A few assumptions must hold for the support domain.
     *   - Each node should be its own first first support node.
     *   - For symmetrical supports the central node may not affect the
     *     calculated derivation and therefore this function is unable to
     *     give a value and will crash.
     *
     * The function is constructed in such a way to make explicit iteration with Neumann
     * boundary conditions simple. To satisfy @f$\dpar{u}{\vec n} = v@f$, we simply use
     * @code
     * u[i] = op.neumann(u, i, domain.normal(i), v);
     * @endcode
     *
     * @param u Scalar field to evaluate.
     * @param node Index of the node, whose value should be changed to satisfy the boundary
     * condition.
     * @param normal Vector @f$\vec{n}@f$ with the same dimensionality as the
     * problem domain giving the direction of the derivation. Usually the unit normal vector.
     * @param val The value that @f$\dpar{ u}{ \vec n}@f$ should be equal to.
     *
     * @return A new value @f$a@f$, such that @f$\dpar{ u}{ \vec n} = v@f$
     * after setting @f$ u(i) = a@f$.
     *
     * @throws Assertion fails if `node`'s first support node is not `node`.
     * @throws Assertion fails if `node`'s value has negligible effect on the derivation, because
     * no value could change the value of @f$\dpar{ u}{ \vec n}@f$.
     */
    template <typename scalar_field_t>
    typename scalar_type<scalar_field_t>::type neumann(
            const scalar_field_t& u, int node, const vector_t& normal,
            typename scalar_type<scalar_field_t>::type val) const;

    /**
     * Returns an approximation of requested derivative of field `u` in `node`-th node.
     * @return The approximation of the derivative @f$ \dpar{ u}{ x_v} @f$,
     * where @f$x_v@f$ is the `var`-th variable and approximation is taken from @ref ss.
     * @throws Assertion fails unless `0 <= var && var < dim` holds.
     */
    template <typename scalar_field_t>
    typename scalar_type<scalar_field_t>::type d1(
            const scalar_field_t& u, int var, int node) const {
        return apply<Der1s<dim>, scalar_field_t>(u, node, Der1<dim>(var));
    }

    /**
     * Returns an approximation of requested second derivative of field `u` in `node`-th node.
     * @return The approximation of the derivative @f$ \dpar{^2 u}{ x_v x_w} @f$,
     * where @f$x_v@f$ is the `varmin`-th variable, @f$x_w@f$ is the `varmax`-th variable
     * and approximation is taken from @ref ss.
     * @throws Assertion fails unless `0 <= varmin <= varmax < dim` holds.
     */
    template <typename scalar_field_t>
    typename scalar_type<scalar_field_t>::type d2(const scalar_field_t& u, int varmin,
                                                  int varmax, int node) const {
        return apply<Der2s<dim>, scalar_field_t>(u, node, Der2<dim>(varmin, varmax));
    }

    /**
     * Returns an approximation of applying the requested operator to field `u` at node `node`.
     * @param u Field to which the operator is applied.
     * @param node Node at which the derived field is evaluated.
     * @param o The operator to apply.
     * @warning The operator approximation is obtained from shape storage using the
     * `op_family_t` type. The first weights computed for family of type `op_family_t` are used.
     *
     * Usage example:
     * @snippet test/operators/ExplicitOperators_test.cpp Custom explicit operators
     */
    template <typename op_family_t, typename scalar_field_t>
    typename scalar_type<scalar_field_t>::type apply(
            const scalar_field_t& u, int node, typename op_family_t::operator_t o) const;


    /// Overload for default constructible operator. @sa apply
    template <typename op_family_t, typename scalar_field_t>
    typename scalar_type<scalar_field_t>::type apply(const scalar_field_t& u, int node) const {
        return apply<op_family_t, scalar_field_t>(u, node, {});
    }

    /// Output basic information about given operators.
    template <typename S>
    friend std::ostream& operator<<(std::ostream& os, const ExplicitOperators<S>& op);
};

}  // namespace mm

#endif  // MEDUSA_BITS_OPERATORS_EXPLICITOPERATORS_FWD_HPP_
