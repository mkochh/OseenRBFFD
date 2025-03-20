#ifndef MEDUSA_BITS_SPATIAL_SEARCH_KDTREE_FWD_HPP_
#define MEDUSA_BITS_SPATIAL_SEARCH_KDTREE_FWD_HPP_

/**
 * @file
 * Declaration of KDTree.
 *
 * @example test/spatial_search/KDTree_test.cpp
 */

#include <medusa/Config.hpp>
#include <nanoflann/nanoflann.hpp>
#include <medusa/bits/types/Range_fwd.hpp>
#include <iosfwd>
#include <medusa/bits/utils/memutils.hpp>
#include "PointCloud.hpp"

namespace mm {

/**
 * Class representing a static [k-d tree](https://en.wikipedia.org/wiki/K-d_tree) data structure.
 *
 * This class is a wrapper around [nanoflann library](https://github.com/jlblancoc/nanoflann)
 * for nearest neighbours. For specific behaviors please refer to
 * [nanoflann documentation](http://jlblancoc.github.io/nanoflann/).
 *
 * If dynamic insertions and removals are needed, use KDTreeMutable.
 *
 * Usage example:
 * @snippet spatial_search/KDTree_test.cpp KDTree usage example
 * @ingroup utils
 *
 * @sa KDTreeMutable
 */
template <class vec_t>
class KDTree {
  public:
    typedef vec_t vector_t;   ///< Vector type used.
    typedef typename vector_t::scalar_t scalar_t;   ///< Scalar type used.
    /// Store dimension of the space.
    enum { /** Dimensionality of the space. */ dim = vec_t::dim };

  private:
    typedef nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<scalar_t, kdtree_internal::PointCloud<vec_t>>,
            kdtree_internal::PointCloud<vec_t>, dim, int> kd_tree_t;   ///< k-d tree type.

    kdtree_internal::PointCloud<vec_t> points_;  ///< Points, contained in the tree.
    kd_tree_t tree;  ///< Actual tree build over points.

  public:
    /**
     * Constructor that builds the search tree for the given points
     * @param points A collection of points.
     */
    explicit KDTree(const Range<vec_t>& points) :
            points_(points), tree(dim, points_, nanoflann::KDTreeSingleIndexAdaptorParams(20)) {
        tree.buildIndex();
    }

    /// Creates an empty KDTree. The tree may later be filled using KDTree::resetTree.
    KDTree() : points_(), tree(dim, points_, nanoflann::KDTreeSingleIndexAdaptorParams(20)) {}

    /**
     * Grows a new tree with new points.
     * This function deletes the old points and the old tree and builds a new one with the given
     * points.
     *
     * @param points A new container of points.
     */
    void reset(const Range<vec_t>& points) {
        points_.setPts(points);
        tree.buildIndex();
    }

    /**
     * Find `k` nearest neighbors to given `point`. This method uses ANN `query` function.
     *
     * @param point Find closest points to this point.
     * @param k How many nearest points to find.
     *
     * @return A pair of two vectors of size `k` containing indices of nearest
     * neighbors and squared distances from `point` to its neighbours.
     */
    std::pair<Range<int>, Range<scalar_t>> query(const vec_t& point, int k = 1) const;

    /**
     * Find neighbors of `point` in given radius.
     *
     * @param point Find closest points to this point.
     * @param radius_squared Maximum distance (squared) to search.
     *
     * @return A pair of two vectors containing indices of nearest
     * neighbors and squared distances from `point` to its neighbours.
     */
    std::pair<Range<int>, Range<scalar_t>> query(
            const vec_t& point, const scalar_t& radius_squared) const;

    /**
     * Get the coordinates of a point in the tree.
     * Given the index of a point as it appeared in the original list
     * this function returns a its coordinates.
     * @note This function is slow as it has to copy the values from its internal
     * containers. It should not be used as a substitute for storing your own
     * values.
     *
     * Example:
     * @snippet KDTree_test.cpp KDTree get
     *
     * @param index Index of the point as given in the constructor.
     * @return Vector object with the coordinates of the `index`-th point.
     */
    vec_t get(int index) const { return points_.get(index); }

    /// Vectorized version of KDTree::get
    Range<vec_t> get(const Range<int>& indexes) const {
        const int n = indexes.size();
        Range<vec_t> result(n);
        for (int i = 0; i < n; ++i) {
            result[i] = points_.get(indexes[i]);
        }
        return result;
    }

    /// Returns number of points in this tree.
    int size() const { return points_.kdtree_get_point_count(); }

    /// Output basic info about given tree.
    template <class V>
    friend std::ostream& operator<<(std::ostream& os, const KDTree<V>& tree);
};

}  // namespace mm

#endif  // MEDUSA_BITS_SPATIAL_SEARCH_KDTREE_FWD_HPP_
