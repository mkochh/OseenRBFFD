#ifndef MEDUSA_BITS_OPERATORS_EXPLICITOPERATORS_HPP_
#define MEDUSA_BITS_OPERATORS_EXPLICITOPERATORS_HPP_

/**
 * @file
 * Explicit operators implementation.
 */

#include "ExplicitOperators_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/types/Vec.hpp>
#include "ShapeStorage_fwd.hpp"

namespace mm {

/**
 * Performs all asserts for a given node: asserts that node index is valid and that shape functions
 * were initialized for a given node.
 * @param node Node index to consider.
 */
#define NODE_ASSERTS(node) \
    assert_msg(0 <= (node) && (node) < ss->size(), "Node index %d must be in range " \
               "[0, num_nodes = %d).", ss->size());

/// @cond
template <typename shape_storage_type>
template <typename op_family_t, typename scalar_field_t>
typename scalar_type<scalar_field_t>::type
ExplicitOperators<shape_storage_type>::apply(const scalar_field_t& u, int node,
                                             typename op_family_t::operator_t o) const {
    NODE_ASSERTS(node);
    int op_idx = op_family_t::index(o);
    typename scalar_type<scalar_field_t>::type ret = 0.0;
    for (int j = 0; j < ss->supportSize(node); ++j)
        ret += ss->template get<op_family_t>(op_idx, node, j) * u[ss->support(node, j)];
    return ret;
}

template <typename shape_storage_type>
template <typename scalar_field_t>
Vec<typename scalar_type<scalar_field_t>::type, ExplicitOperators<shape_storage_type>::dim>
ExplicitOperators<shape_storage_type>::grad(const scalar_field_t& u, int node) const {
    NODE_ASSERTS(node);
    Vec<typename scalar_type<scalar_field_t>::type, ExplicitOperators<shape_storage_type>::dim> ret;
    for (int var = 0; var < dim; ++var) {
        ret[var] = ss->d1(var, node, 0) * u[ss->support(node, 0)];
        for (int i = 1; i < ss->supportSize(node); ++i) {
            ret[var] += ss->d1(var, node, i) * u[ss->support(node, i)];
        }
    }
    return ret;
}


template <typename shape_storage_type>
template <typename scalar_field_t>
typename scalar_type<scalar_field_t>::type
ExplicitOperators<shape_storage_type>::neumann(
        const scalar_field_t& u, int node, const vector_t& normal,
        typename scalar_type<scalar_field_t>::type val) const {
    NODE_ASSERTS(node);
    assert_msg(ss->support(node, 0) == node, "First support node should be the node itself.");
    scalar_t denominator = 0;
    for (int d = 0; d < dim; ++d) {
        for (int i = 1; i < ss->supportSize(node); ++i) {
            val -= normal[d] * ss->d1(d, node, i) * u[ss->support(node, i)];
        }
        // i = 0
        denominator += normal[d] * ss->d1(d, node, 0);
    }
    assert_msg(std::abs(denominator) >= 1e-14,
               "Node %d has no effect on the flux in direction %s. The cause of this might be wrong"
               " normal direction, bad neighbourhood choice or bad nodal positions.", node, normal);
    return val / denominator;
}
/// @endcond

/// Output basic information about given operators.
template <typename S>
std::ostream& operator<<(std::ostream& os, const ExplicitOperators<S>& op)  {
    if (!op.hasShapes()) {
        return os << "Explicit operators without any linked storage.";
    }
    return os << "Explicit operators over storage: " << *op.ss;
}

template <typename Derived, typename vec_t, typename OpFamilies>
ExplicitOperators<Derived>
ShapeStorage<Derived, vec_t, OpFamilies>::explicitOperators() const {
    return ExplicitOperators<Derived>(*static_cast<const Derived*>(this));
}

}  // namespace mm

#endif  // MEDUSA_BITS_OPERATORS_EXPLICITOPERATORS_HPP_
