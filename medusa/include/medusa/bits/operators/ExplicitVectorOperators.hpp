#ifndef MEDUSA_BITS_OPERATORS_EXPLICITVECTOROPERATORS_HPP_
#define MEDUSA_BITS_OPERATORS_EXPLICITVECTOROPERATORS_HPP_

/**
 * @file
 * Explicit vector operators implementation.
 */

#include "ExplicitVectorOperators_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <Eigen/Core>
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
template <class shape_storage_type>
template <class vector_field_t>
Eigen::Matrix<
        typename vector_type<vector_field_t>::type::scalar_t,
        vector_type<vector_field_t>::type::dim,
        ExplicitVectorOperators<shape_storage_type>::dim>
ExplicitVectorOperators<shape_storage_type>::grad(const vector_field_t& u, int node) const {
    NODE_ASSERTS(node);
    constexpr int vector_dim = vector_type<vector_field_t>::type::dim;
    constexpr int domain_dim = ExplicitVectorOperators<shape_storage_type>::dim;
    Eigen::Matrix<typename vector_type<vector_field_t>::type::scalar_t, vector_dim, domain_dim> ret;
    for (int d = 0; d < vector_dim; ++d) {
        for (int var = 0; var < domain_dim; ++var) {
            ret(d, var) = ss->d1(var, node, 0) * u[ss->support(node, 0)][d];
            for (int i = 1; i < ss->supportSize(node); ++i) {
                ret(d, var) += ss->d1(var, node, i) * u[ss->support(node, i)][d];
            }
        }
    }
    return ret;
}

template <class shape_storage_type>
template <typename op_family_t, typename vector_field_t>
typename vector_type<vector_field_t>::type
ExplicitVectorOperators<shape_storage_type>::apply(const vector_field_t& u, int node,
                                                   typename op_family_t::operator_t o) const {
    constexpr int u_dim = vector_type<vector_field_t>::type::dim;
    typename vector_type<vector_field_t>::type res;
    res.setZero();
    int idx = op_family_t::index(o);
    for (int i = 0; i < ss->supportSize(node); ++i) {
        for (int d = 0; d < u_dim; ++d) {
            res[d] += ss->template get<op_family_t>(idx, node, i) * u[ss->support(node, i)][d];
        }
    }
    return res;
}

template <class shape_storage_type>
template <class vector_field_t>
typename vector_type<vector_field_t>::type::scalar_t
ExplicitVectorOperators<shape_storage_type>::div(const vector_field_t& u, int node) const {
    NODE_ASSERTS(node);
    static_assert(static_cast<int>(dim) == static_cast<int>(vector_type<vector_field_t>::type::dim),
                  "Domain and filed dimensions must match.");
    typename vector_type<vector_field_t>::type::scalar_t ret =
            ss->d1(0, node, 0) * u[ss->support(node, 0)][0];
    for (int i = 1; i < ss->supportSize(node); ++i)
        ret += ss->d1(0, node, i) * u[ss->support(node, i)][0];
    for (int var = 1; var < dim; ++var)
        for (int i = 0; i < ss->supportSize(node); ++i)
            ret += ss->d1(var, node, i) * u[ss->support(node, i)][var];
    return ret;
}

template <class shape_storage_type>
template <class vector_field_t>
typename vector_type<vector_field_t>::type
ExplicitVectorOperators<shape_storage_type>::curl(const vector_field_t& u, int node) const {
    NODE_ASSERTS(node);
    static_assert(static_cast<int>(vector_type<vector_field_t>::type::dim) == 3,
            "The field must be three dimensional.");
    typename vector_type<vector_field_t>::type ret;
    ret.setZero();
    for (int i = 0; i < ss->supportSize(node); ++i) {
        ret[0] += ss->d1(1, node, i) * u[ss->support(node, i)][2]
                - ss->d1(2, node, i) * u[ss->support(node, i)][1];
        ret[1] += ss->d1(2, node, i) * u[ss->support(node, i)][0]
                - ss->d1(0, node, i) * u[ss->support(node, i)][2];
        ret[2] += ss->d1(0, node, i) * u[ss->support(node, i)][1]
                - ss->d1(1, node, i) * u[ss->support(node, i)][0];
    }
    return ret;
}

template <class shape_storage_type>
template <class vector_field_t>
typename vector_type<vector_field_t>::type
ExplicitVectorOperators<shape_storage_type>::graddiv(const vector_field_t& u, int node) const {
    NODE_ASSERTS(node);
    static_assert(static_cast<int>(dim) == static_cast<int>(vector_type<vector_field_t>::type::dim),
                  "Domain and filed dimensions must match.");
    typename vector_type<vector_field_t>::type ret;
    ret.setZero();
    for (int d1 = 0; d1 < dim; ++d1) {
        for (int d2 = 0; d2 < dim; ++d2) {  // loop over dimensions
            int dmin = std::min(d1, d2);
            int dmax = std::max(d1, d2);
            for (int i = 0; i < ss->supportSize(node); ++i) {
                ret[d1] += ss->d2(dmin, dmax, node, i) * u[ss->support(node, i)][d2];
            }
        }
    }
    return ret;
}

template <class shape_storage_type>
template <class vector_field_t>
typename vector_type<vector_field_t>::type
ExplicitVectorOperators<shape_storage_type>::neumann(
        const vector_field_t& u, int node, const vector_t& normal,
        typename vector_type<vector_field_t>::type val) const {
    NODE_ASSERTS(node);
    static_assert(static_cast<int>(dim) == static_cast<int>(vector_type<vector_field_t>::type::dim),
            "Domain and filed dimensions must match.");
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

/// Output basic info about given operators.
template <typename S>
std::ostream& operator<<(std::ostream& os, const ExplicitVectorOperators<S>& op) {
    if (!op.hasShapes()) {
        return os << "Explicit vector operators without any linked storage.";
    }
    return os << "Explicit vector operators over storage: " << *op.ss;
}

template <typename Derived, typename vec_t, typename OpFamilies>
ExplicitVectorOperators<Derived>
ShapeStorage<Derived, vec_t, OpFamilies>::explicitVectorOperators() const {
    return ExplicitVectorOperators<Derived>(*static_cast<const Derived*>(this));
}

}  // namespace mm

#endif  // MEDUSA_BITS_OPERATORS_EXPLICITVECTOROPERATORS_HPP_
