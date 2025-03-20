#ifndef MEDUSA_BITS_SPATIAL_SEARCH_POINTCLOUD_HPP_
#define MEDUSA_BITS_SPATIAL_SEARCH_POINTCLOUD_HPP_

/**
 * @file
 * Implementation of KDTree storage class.
 */

#include <vector>

namespace mm {

/// Implementation details of KDTree, contains intermediate storage class.
namespace kdtree_internal {

/// Helper class for KDTree with appropriate accessors containing a set of points. @ingroup utils
template <typename vec_t>
struct PointCloud {
    std::vector<vec_t> pts;  ///< Points, contained in the tree.

    /// Construct an empty point set.
    PointCloud() = default;

    /// Construct from an array of points.
    PointCloud(const std::vector<vec_t>& pts) : pts(pts) {}

    /// Reset contained points.
    void setPts(const std::vector<vec_t>& pts) {
        PointCloud::pts = pts;
    }

    /// Interface requirement: returns number of data points.
    inline int kdtree_get_point_count() const { return pts.size(); }

    /// Interface requirement: returns `dim`-th coordinate of `idx`-th point.
    inline typename vec_t::scalar_t kdtree_get_pt(const size_t idx, int dim) const {
        return pts[idx][dim];
    }

    /// Access the points.
    inline vec_t get(const size_t idx) const { return pts[idx]; }

    /// Add a point to the cloud.
    inline void add(const vec_t& p) { pts.push_back(p); }

    /// Comply with the interface.
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};

}  // namespace kdtree_internal
}  // namespace mm

#endif  // MEDUSA_BITS_SPATIAL_SEARCH_POINTCLOUD_HPP_
