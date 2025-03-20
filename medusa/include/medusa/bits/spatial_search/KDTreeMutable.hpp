#ifndef MEDUSA_BITS_SPATIAL_SEARCH_KDTREEMUTABLE_HPP_
#define MEDUSA_BITS_SPATIAL_SEARCH_KDTREEMUTABLE_HPP_

/**
 * @file
 * Implementation of KDTreeMutable.
 */

#include "KDTreeMutable_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <cmath>

namespace mm {

/// Output basic info about given tree.
template <class V>
std::ostream& operator<<(std::ostream& os, const KDTreeMutable<V>& tree)  {
    return os << "KDTreeMutable:\n"
              << "    dimension: " << tree.dim << '\n'
              << "    num of points: " << tree.size() << std::endl;
}

template <class vec_t>
void KDTreeMutable<vec_t>::insert(const vec_t& point) {
    assert_msg(point.array().isFinite().prod() == 1, "Invalid point.");
    auto n = points_.kdtree_get_point_count();
    points_.add(point);
    tree.addPoints(n, n);
    ++size_;
}

template <class vec_t>
void KDTreeMutable<vec_t>::insert(const Range<vec_t>& points) {
    auto n = points_.kdtree_get_point_count();
    for (const auto& p : points) {
        assert_msg(p.array().isFinite().prod() == 1, "One of the points is invalid.");
        points_.add(p);
    }
    size_ += points.size();
    tree.addPoints(n, n + points.size() - 1);
}

template <class vec_t>
std::pair<mm::Range<int>, mm::Range<double>>
KDTreeMutable<vec_t>::query(const vec_t& point, int k) {
    assert_msg(point.array().isFinite().prod() == 1, "Invalid query point.");
    nanoflann::KNNResultSet<scalar_t, int> resultSet(k);
    Range<int> ret_index(k);
    Range<scalar_t> out_dist_sqr(k);
    resultSet.init(&ret_index[0], &out_dist_sqr[0]);
    tree.findNeighbors(resultSet, point.data(), nanoflann::SearchParams(k));
    assert_msg(resultSet.full(), "Not enough points in the tree, you requested %d points, "
                                 "but the tree contains only %d points.", k, size());
    return {ret_index, out_dist_sqr};
}

}  // namespace mm

#endif  // MEDUSA_BITS_SPATIAL_SEARCH_KDTREEMUTABLE_HPP_
