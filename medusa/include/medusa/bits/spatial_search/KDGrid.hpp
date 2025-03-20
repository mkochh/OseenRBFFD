#ifndef MEDUSA_BITS_SPATIAL_SEARCH_KDGRID_HPP_
#define MEDUSA_BITS_SPATIAL_SEARCH_KDGRID_HPP_

/**
 * @file
 * Implementation of KDGrid.
 */

#include "KDGrid_fwd.hpp"
#include "Grid.hpp"
#include <medusa/bits/utils/numutils.hpp>
#include <medusa/bits/utils/assert.hpp>
#include <iostream>

namespace mm {

template <typename vec_t>
typename KDGrid<vec_t>::IndexArray
KDGrid<vec_t>::compute_size(const vec_t& bot, const vec_t& top, scalar_t cell_size) {
    IndexArray counts;
    for (int i = 0; i < dim; ++i) {
        assert_msg(bot[i] <= top[i], "Bottom bounds %s must be below top bounds %s.", bot, top);
        counts[i] = iceil((top[i] - bot[i]) / cell_size) + 1;
    }
    return counts;
}

template <typename vec_t>
bool KDGrid<vec_t>::increment(KDGrid::IndexArray& cur, int span, const KDGrid::IndexArray& ref) {
    for (int i = 0; i < dim; ++i) {
        if (cur[i] == ref[i] + span) {
            cur[i] = ref[i] - span;
        } else {
            cur[i] += 1;
            return true;
        }
    }
    return false;
}

template <typename vec_t>
int KDGrid<vec_t>::insert(const vec_t& p, const KDGrid::IndexArray& index) {
    assert_msg(grid_(index) == -1, "This grid cell is already occupied by point %d. Consider"
                                   " using smaller cell size.", grid_(index));
    int idx = points_.size();
    grid_(index) = idx;
    points_.push_back(p);
    return idx;
}

template <typename vec_t>
bool
KDGrid<vec_t>::existsPointInSphere(const vec_t& p, scalar_t r, const KDGrid::IndexArray& index) {
    int span = static_cast<int>(std::ceil(r/cell_size_));
    scalar_t r2 = r*r;
    IndexArray tmp_index = index;
    for (int i = 0; i < dim; ++i) tmp_index[i] -= span;
    do {
        if (grid_.inBounds(tmp_index) && grid_(tmp_index) != -1) {
            auto gp = points_[grid_(tmp_index)];
            if ((p-gp).squaredNorm() < r2) {
                return true;
            }
        }
    } while (increment(tmp_index, span, index));
    return false;
}

template <typename vec_t>
typename KDGrid<vec_t>::IndexArray KDGrid<vec_t>::compute_index(const vec_t& p) {
    Eigen::Matrix<scalar_t, dim, 1> idx = (p-bot_) / cell_size_;
    IndexArray index;
    for (int i = 0; i < dim; ++i) {
        index[i] = static_cast<int>(idx[i]);
    }
    return index;
}

/// @cond
template <typename V>
std::ostream& operator<<(std::ostream& os, const KDGrid<V>& search) {
    os << "KDGrid search structure from " << search.bot_ << " to " << search.top_
       << " with cell_size " << search.cell_size_ << " containing " << search.size()
       << " points";
    return os;
}
/// @endcond

}  // namespace mm

#endif  // MEDUSA_BITS_SPATIAL_SEARCH_KDGRID_HPP_
