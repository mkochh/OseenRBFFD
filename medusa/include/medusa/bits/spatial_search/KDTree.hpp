#ifndef MEDUSA_BITS_SPATIAL_SEARCH_KDTREE_HPP_
#define MEDUSA_BITS_SPATIAL_SEARCH_KDTREE_HPP_

/**
 * @file
 * Implementation of KDTree.
 */

#include "KDTree_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/utils/memutils.hpp>

namespace mm {

template <class vec_t>
std::pair<Range<int>, Range<typename vec_t::scalar_t>> KDTree<vec_t>::query(const vec_t& point,
                                                                            int k) const {
    assert_msg(point.array().isFinite().prod() == 1, "Invalid point.");
    Range<int> ret_index(k);
    Range<scalar_t> out_dist_sqr(k);
    int actual_k = tree.knnSearch(point.data(), k, &ret_index[0], &out_dist_sqr[0]);
    assert_msg(actual_k == k, "There were not enough points in the tree, you requested %d "
                              "points, the tree only contains %d points.", k, actual_k);
    return {ret_index, out_dist_sqr};
}

template <class vec_t>
std::pair<Range<int>, Range<typename vec_t::scalar_t>> KDTree<vec_t>::query(
        const vec_t& point, const scalar_t& radius_squared) const {
    assert_msg(point.array().isFinite().prod() == 1, "Invalid point.");
    std::vector<std::pair<int, scalar_t>> idx_dist;
    int k = tree.radiusSearch(point.data(), radius_squared, idx_dist, nanoflann::SearchParams());
    Range<int> idx(k); Range<scalar_t> dists(k);
    for (int i = 0; i < k; ++i) {
        std::tie(idx[i], dists[i]) = idx_dist[i];
    }
    return {idx, dists};
}

/// Output basic info about given tree.
template <class V>
std::ostream& operator<<(std::ostream& os, const KDTree<V>& tree) {
    return os << "KDTree:\n"
              << "    dimension: " << tree.dim << '\n'
              << "    num of points: " << tree.size() << "\n"
              << "    first point: " << tree.get(0) << "\n"
              << "    last point: " << tree.get(tree.size() - 1) << "\n";
}

}  // namespace mm

#endif  // MEDUSA_BITS_SPATIAL_SEARCH_KDTREE_HPP_
