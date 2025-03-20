#ifndef MEDUSA_BITS_OPERATORS_EXPLICITVECTOROPERATORS_FWD_HPP_
#define MEDUSA_BITS_OPERATORS_EXPLICITVECTOROPERATORS_FWD_HPP_

/**
 * @file
 * Explicit vector operators declarations.
 *
 * @example test/operators/ExplicitVectorOperators_test.cpp
 */

#include <medusa/Config.hpp>
#include <Eigen/Core>
#include <medusa/bits/types/traits.hpp>
#include <medusa/bits/approximations/Operators_fwd.hpp>

namespace mm {

/**
 * A class for evaluating typical operators needed in spatial discretization.
 * This class represents vector differential operators, discretely approximated
 * with shape functions that are stored in given shape storage. These shape functions
 * are used to construct approximations of common operators such as gradient, Laplacian,
 * divergence, ... that can be applied to vector fields at given points.
 *
 * @note The scalar types (and sometimes even dimensions) of the given fields and
 * the shapes can be different. We support computing the gradient of a `3d` complex-valued
 * field over 2d domain with `double` shapes: the result is a `3x2` complex matrix.
 *
 * @tparam shape_storage_type Storage class for shapes. Must satisfy the @ref ss-concept.
 *
 * Usage example:
 * @snippet operators/ExplicitVectorOperators_test.cpp Explicit vector operators usage example
 * @ingroup operators
 *
 * @sa ExplicitOperators
 */
template <class shape_storage_type>
class ExplicitVectorOperators {
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
    /// Construct empty explicit vector operators.
    ExplicitVectorOperators() : ss(nullptr) {}
    /**
     * Construct explicit vector operators over given shape storage.
     * The pointer to given storage is stored as a non-owning pointer. It is the user's
     * responsibility that given storage is valid throughout the use use of this class.
     *
     * This class is usually constructed directly from shape storage using the
     * `explicitVectorOperators()` member function.
     */
    explicit ExplicitVectorOperators(const shape_storage_t& ss) : ss(&ss) {}

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
     * Returns Laplacian of vector field `u` in node with index `node`.
     * @returns The Laplacian of `u`, e.g.\ when `dim = 3` and @f$\vec{u} = [u, v, w]@f$, the
     * returned value is
     * @f$ \nabla^2 \vec{u} = [
     * \dpar{^2 u}{x^2} + \dpar{^2 u}{y^2} + \dpar{^2 u}{z^2},
     * \dpar{^2 v}{x^2} + \dpar{^2 v}{y^2} + \dpar{^2 v}{z^2},
     * \dpar{^2 w}{x^2} + \dpar{^2 w}{y^2} + \dpar{^2 w}{z^2}
     * ]^T@f$ evaluated in point with index `node`
     * where derivatives are approximated using shape functions from @ref ss.
     * @throws Static-assertion fails if the dimension of `u` is not equal to domain dimension.
     */
    template <class vector_field_t>
    typename vector_type<vector_field_t>::type lap(const vector_field_t& u, int node) const {
        return apply<Lap<dim>, vector_field_t>(u, node, {});
    }

    /**
     * Returns an approximation of applying the requested operator to each component of
     * field `u` at node `node`.
     * @param u Field to which the operator is applied.
     * @param node Node at which the derived field is evaluated.
     * @param o The operator to apply.
     * @warning The operator approximation is obtained from shape storage using the
     * `op_family_t` type. The first weights computed for family of type `op_family_t` are used.
     *
     * Usage example:
     * @snippet test/operators/ExplicitVectorOperators_test.cpp Custom vector operators
     */
    template <typename op_family_t, typename vector_field_t>
    typename vector_type<vector_field_t>::type apply(
            const vector_field_t& u, int node, typename op_family_t::operator_t o) const;


    /// Overload for default constructible operator. @sa apply
    template <typename op_family_t, typename vector_field_t>
    typename vector_type<vector_field_t>::type apply(const vector_field_t& u, int node) const {
        return apply<op_family_t, vector_field_t>(u, node, {});
    }


    /**
     * Returns gradient of vector field `u` in node with index `node`. The field `u` and
     * the domain can have different dimensions.
     * @returns The gradient (Jacobi matrix) of `u`, e.g.\ when `dim(u) = 3` and
     * and @f$\vec{u} = [u, v, w]@f$ and the domain dimension is 2, the
     * returned value is
     * @f[ \nabla \vec{u} = \begin{bmatrix}
     * \dpar{u}{x}, \dpar{u}{y} \\
     * \dpar{v}{x}, \dpar{v}{y} \\
     * \dpar{w}{x}, \dpar{w}{y}
     * \end{bmatrix} @f]  evaluated in point with index `node`,
     * where derivatives are approximated using shape functions from @ref ss.
     */
    template <class vector_field_t>
    Eigen::Matrix<typename vector_type<vector_field_t>::type::scalar_t,
                  vector_type<vector_field_t>::type::dim, dim>
    grad(const vector_field_t& u, int node) const;

    /**
     * Returns divergence of a vector field `u` in in node with index `node`.
     * @returns The divergence of `u`, e.g.\ when `dim = 3` and @f$\vec{u} = [u, v, w]@f$, the
     * returned value is
     * @f$ \nabla\cdot \vec{u} = \dpar{u}{x} + \dpar{v}{y} + \dpar{w}{z} @f$
     * evaluated in point with index `node`,
     * where derivatives are approximated using shape functions from @ref ss.
     * @throws Static-assertion fails if the dimension of `u` is not equal to domain dimension.
     */
    template <class vector_field_t>
    typename vector_type<vector_field_t>::type::scalar_t div(
            const vector_field_t& u, int node) const;

    /**
     * Returns curl of a vector field `u` in node with index `node`.
     * @return The rotor of `u`, where `dim = 3` and @f$\vec{u} = [u, v, w]@f$. The return value is
     * @f$ \curl \vec{u} = \[\dpar{w}{y} - \dpar{v}{z}, \dpar{u}{z} - \dpar{w}{x}, \dpar{v}{x} - \dpar{u}{y}\] @f$
     * evaluated in point with index `node`, where derivatives are approximated using shape functions from @ref ss.
     * @throws Static-assertion fails if `u` is not three dimensional.
     */
     template <class vector_field_t>
     typename vector_type<vector_field_t>::type curl(const vector_field_t& u, int node) const;

    /**
     * Returns gradient of divergence of a vector field `u` in in node with index `node`.
     * @returns The gradient of divergence of `u`, e.g.\ when `dim = 3` and
     * @f$\vec{u} = [u, v, w]@f$, the returned value is
     * @f[ \nabla(\nabla\cdot \vec{u}) = \begin{bmatrix}
     *   \dpar{^2u}{x^2} + \dpar{^2v}{x\partial y} + \dpar{^2w}{x\partial z} \\
     *   \dpar{^2u}{x\partial y} + \dpar{^2v}{^2y} + \dpar{^2w}{y\partial z} \\
     *   \dpar{^2u}{x\partial z} + \dpar{^2v}{y\partial z} + \dpar{^2w}{^2z}
     * \end{bmatrix} @f]
     * evaluated in point with index `node`,
     * where derivatives are approximated using shape functions from @ref ss.
     * @throws Static-assertion fails if the dimension of `u` is not equal to domain dimension.
     */
    template <class vector_field_t>
    typename vector_type<vector_field_t>::type graddiv(const vector_field_t& u, int node) const;
    /**
     * This is a vectorised version of ExplicitOperators::neumann.
     * Calculates a new field value `a` for `u` at `node`, such that the Neumann boundary
     * condition `du/dn(node) = val` holds.
     * Denote given index `node` by @f$i@f$ and `val` by @f$\vec v@f$.
     * The function returns such @f$\vec a@f$, that by setting @f$ \vec{u}(i) = \vec{a}@f$ we
     * satisfy the condition @f$ \dpar{\vec u}{\vec n}(i) = \vec v@f$. This is done using the
     * following formula:
     * @f[ \vec{a} = \frac{ \vec{v} -
     *    \sum\limits_{i=1}^{s} \sum\limits_{d=1}^{dim}
     *      n_d \chi_{\dpar{}{ x_d}}(i)\ \vec{u}(i)
     * }{\sum\limits_{d=1}^{dim} n_d \chi_{\dpar{}{ x_d}}(0)},
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
     * @param u Scalar field to evaluate.
     * @param node Index of the node, whose value should be changed to satisfy the boundary
     * condition.
     * @param normal Vector @f$\vec{n}@f$ with the same dimensionality as the
     * problem domain giving the direction of the derivation. Usually the unit normal vector.
     * @param val The value that @f$\dpar{ \vec u}{ \vec n}@f$ should be equal to.
     *
     * @return A new value @f$\vec a@f$, such that @f$\dpar{\vec u}{\vec n} = \vec v@f$
     * after setting @f$\vec u(i) = \vec a@f$.
     *
     * @throws Assertion fails if `node`'s first support node is not `node`.
     * @throws Assertion fails if `node`'s value has negligible effect on the derivation, because
     * no value could change the value of @f$\dpar{\vec u}{\vec n}@f$.
     * @throws Static-assertion fails if the dimension of `u` is not equal to domain dimension.
     *
     * @sa ExplicitOperators::neumann
     */
    template <class vector_field_t>
    typename vector_type<vector_field_t>::type neumann(const vector_field_t& u, int node,
            const vector_t& normal, typename vector_type<vector_field_t>::type val) const;

    /// Output basic info about given operators.
    template <typename S>
    friend std::ostream& operator<<(std::ostream& os, const ExplicitVectorOperators<S>& op);
};

}  // namespace mm

#endif  // MEDUSA_BITS_OPERATORS_EXPLICITVECTOROPERATORS_FWD_HPP_
