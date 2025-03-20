#ifndef MEDUSA_BITS_OPERATORS_COMPUTESHAPES_HPP_
#define MEDUSA_BITS_OPERATORS_COMPUTESHAPES_HPP_

/**
 * @file
 * Implementations of shape function computation utilities.
 */

#include "computeShapes_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/domains/DomainDiscretization.hpp>
#include <medusa/bits/utils/stdtypesutils.hpp>
#include <Eigen/Core>
#include <tuple>
#include <type_traits>

namespace mm {

/// Shortcut to compute shapes for given domain from operator mask.
template <class vec_t>
template <sh::shape_flags mask, typename approx_t>
RaggedShapeStorage<vec_t, typename sh::operator_tuple<mask, vec_t::dim>::type>
DomainDiscretization<vec_t>::computeShapes(approx_t approx, const indexes_t& indexes) const {
    static_assert(static_cast<int>(vec_t::dim) == static_cast<int>(approx_t::dim),
                  "Domain and approximation engine dimensions do not match");
    RaggedShapeStorage<vec_t, typename sh::operator_tuple<mask, vec_t::dim>::type> storage;
    storage.resize(supportSizes());
    typename sh::operator_tuple<mask, vec_t::dim>::type operators{};
    if (indexes.empty()) mm::computeShapes(*this, approx, all(), operators , &storage);
    else mm::computeShapes(*this, approx, indexes, operators , &storage);
    return storage;
}

/// Shortcut to compute shapes for given domain from operator tuple.
template <class vec_t>
template <typename OperatorTuple, typename approx_t>
RaggedShapeStorage<vec_t, OperatorTuple>
DomainDiscretization<vec_t>::computeShapes(approx_t approx, const indexes_t& indexes) const {
    static_assert(static_cast<int>(vec_t::dim) == static_cast<int>(approx_t::dim),
                  "Domain and approximation engine dimensions do not match");
    RaggedShapeStorage<vec_t, OperatorTuple> storage;
    OperatorTuple operators{};
    storage.resize(supportSizes());
    if (indexes.empty()) mm::computeShapes(*this, approx, all(), operators , &storage);
    else mm::computeShapes(*this, approx, indexes, operators , &storage);
    return storage;
}

/// Shortcut to compute shapes for given domain from given operators.
template <class vec_t>
template <typename OperatorTuple, typename approx_t>
RaggedShapeStorage<vec_t, OperatorTuple>
DomainDiscretization<vec_t>::computeShapes(OperatorTuple operators, approx_t approx,
const indexes_t& indexes) const {
    static_assert(static_cast<int>(vec_t::dim) == static_cast<int>(approx_t::dim),
                  "Domain and approximation engine dimensions do not match");
    RaggedShapeStorage<vec_t, OperatorTuple> storage;
    storage.resize(supportSizes());
    if (indexes.empty()) mm::computeShapes(*this, approx, all(), operators , &storage);
    else mm::computeShapes(*this, approx, indexes, operators , &storage);
    return storage;
}

namespace shapes_internal {

/// @cond
template <typename approx_t, typename shape_storage_t, typename OperatorTuple, std::size_t N>
struct ShapeComputer {
    static void compute(const OperatorTuple& operators, const approx_t& approx, int node,
                        shape_storage_t* storage) {
        // do computation for previous ones
        ShapeComputer<approx_t, shape_storage_t, OperatorTuple, N-1>::compute(
                operators, approx, node, storage);

        // Do computation for current operator family:
        // get family type (for storage indexing)
        typedef typename std::tuple_element<N-1, OperatorTuple>::type op_family_t;
        // get family value, to produce operators
        auto op_family = std::get<N-1>(operators);
        // count operators (requirements ensure compatibility with `index` method used in `apply`)
        int i = 0;
        for (const auto op : op_family.operators()) {
            // Access storage by family type, not by family index (e.g. using N-1), since given
            // operators may only be a subset of total operators in storage, and using
            // N-1 is wrong.
            storage->template setShape<op_family_t>(i, node, approx.getShape(op));
            ++i;
        }
    }
};

template <typename approx_t, typename shape_storage_t, typename OperatorTuple>
struct ShapeComputer<approx_t, shape_storage_t, OperatorTuple, 0> {
    static void compute(const OperatorTuple&, const approx_t&, int, shape_storage_t*) {}
};

}  // namespace shapes_internal
/// @endcond

template <typename approx_t,  typename shape_storage_t, typename ...Operators>
void computeShapes(const DomainDiscretization<typename approx_t::vector_t>& domain,
                   approx_t approx, const indexes_t& indexes,
                   const std::tuple<Operators...>& operators, shape_storage_t* storage) {
    int N = domain.size();
    assert_msg(domain.supports().size() == N,
               "domain.support.size = %d and domain.size = %d, but should be the same. "
               "Did you forget to find support before computing shapes?",
               domain.supports().size(), domain.size());
    for (const auto& c : indexes) {
        assert_msg(0 <= c && c < N,
                   "Index %d is not a valid index of a point in the domain, must be in range "
                   "[%d, %d).", c, 0, N);
        assert_msg(!domain.support(c).empty(),
                   "Node %d has empty support! Did you forget to find support before "
                   "computing shapes?", c);
    }
    assert_msg(storage->size() == N,
               "Storage must of appropriate size before calling computeShapes(), got size %d, "
               "expected %d.", storage->size(), N);

    // construct shape functions for every point specified
    int cc;
    int isize = static_cast<int>(indexes.size());
    // create local copies of mls for each thread
    #if !defined(_OPENMP)
    approx_t& local_approx = approx;
    #else
    std::vector<approx_t> approx_copies(omp_get_max_threads(), approx);
    #pragma omp parallel for private(cc) schedule(static)
    #endif
    for (cc = 0; cc < isize; ++cc) {
        int node = indexes[cc];
        //  store local copy of approximation engine -- for parallel computing
        #if defined(_OPENMP)
        approx_t& local_approx = approx_copies[omp_get_thread_num()];
        #endif
        assert_msg(storage->supportSize(node) == domain.supportSize(node),
                   "Storage must of appropriate size before calling computeShapes(). "
                   "Got support size %d at node %d, expected %d.",
                   storage->supportSize(node), node, domain.supportSize(node));
        // Preps local 1D support domain vector (cache friendly storage).
        storage->setSupport(node, domain.support(node));
        // Resets local approximation engine with local parameters.
        Range<typename approx_t::vector_t> supp_domain = domain.supportNodes(node);
        local_approx.compute(domain.pos(node), supp_domain);
        // Iterate over all operator families to compute approximations.
        shapes_internal::ShapeComputer<
                approx_t, shape_storage_t, std::tuple<Operators...>, sizeof...(Operators)
        >::compute(operators, local_approx, node, storage);
    }
}

}  // namespace mm

#endif  // MEDUSA_BITS_OPERATORS_COMPUTESHAPES_HPP_
