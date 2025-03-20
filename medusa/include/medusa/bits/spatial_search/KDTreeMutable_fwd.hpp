#ifndef MEDUSA_BITS_SPATIAL_SEARCH_KDTREEMUTABLE_FWD_HPP_
#define MEDUSA_BITS_SPATIAL_SEARCH_KDTREEMUTABLE_FWD_HPP_

/**
 * @file
 * Declaration of KDTreeMutable.
 *
 * @example test/spatial_search/KDTreeMutable_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/types/Range_fwd.hpp>
#include <nanoflann/nanoflann.hpp>
#include <array>
#include <iosfwd>
#include "PointCloud.hpp"

namespace mm {

/**
 * A [k-d tree](https://en.wikipedia.org/wiki/K-d_tree) data structure that supports dynamic
 * insertions and lazy-removal.
 *
 * This class is a wrapper around [nanoflann](https://github.com/jlblancoc/nanoflann)
 * k-d tree implementation. For specific behavior please refer to
 * [nanoflann docs](http://jlblancoc.github.io/nanoflann/classnanoflann_1_1KDTreeSingleIndexDynamicAdaptor.html).
 *
 * If dynamic structure is not needed, a faster static KDTree structure can be used.
 *
 * Usage example:
 * @snippet spatial_search/KDTreeMutable_test.cpp KDTreeMutable usage example
 * @ingroup utils
 * @sa KDTree, KDGrid, GeneralFill, FindSupport
 */
template <class vec_t>
class KDTreeMutable {
  public:
    typedef vec_t vector_t;   ///< Vector type used.
    typedef typename vector_t::scalar_t scalar_t;   ///< Scalar type used.
    /// Store dimension of the space.
    enum { /** Dimensionality of the space. */ dim = vec_t::dim };
    static_assert(std::is_same<double, scalar_t>::value,
                  "Not Implemented Error: only implemented for double type!");
  private:
    typedef nanoflann::KDTreeSingleIndexDynamicAdaptor<
            nanoflann::L2_Simple_Adaptor<scalar_t, kdtree_internal::PointCloud<vec_t>>,
            kdtree_internal::PointCloud<vec_t>, dim> kd_tree_t;  ///< Tree type.

    kdtree_internal::PointCloud<vec_t> points_;  ///< Points, contained in the tree.
    int size_;  ///< Number of points in the tree.
    kd_tree_t tree;  ///< Actual tree build over points.

  public:
    /**
    * Constructor that builds the search tree for the given points
    * @param points A collection of points.
    */
    explicit KDTreeMutable(const Range<vec_t>& points) :
            points_(points), size_(points.size()),
            tree(dim, points_, nanoflann::KDTreeSingleIndexAdaptorParams(20)) {}

    /// Creates an empty k-d tree.
    KDTreeMutable() : points_(), size_(0),
            tree(dim, points_, nanoflann::KDTreeSingleIndexAdaptorParams(20)) {}

    /**
     * Grows a new tree with new points.
     * This function deletes the old points and the old tree and builds a new one with the given
     * points.
     *
     * @param points A new container of points.
     */
    void reset(const Range<vec_t>& points) {
        points_.setPts(points);
        size_ = points.size();
        tree.reset();
    }
    /**
     * Insert a point into the tree.
     * @param point Point to be inserted into the tree.
     */
    void insert(const vec_t& point);

    /// Inserts a sequence of points.
    void insert(const Range<vec_t>& points);

    /// Check if any point exists in sphere centered at `p` with radius `r`.
    bool existsPointInSphere(const vec_t& p, scalar_t r) {
        if (size_ == 0) return false;
        return query(p).second[0] <= r*r;
    }

    /**
     * Removes a point with given index from the tree. The indexes of the points are given
     * sequentially at insertion and do not change. The removal is lazy and point still takes
     * up memory.
     */
    void remove(int index) {
        size_ -= tree.removePoint(index);
    }

    /// Returns number of points in the tree.
    int size() const { return size_; }

    /**
     * Find `k` nearest neighbors to given point.
     * @param point Find closest points to this point.
     * @param k How many nearest points to find.
     *
     * @return A pair of two vectors of size `k` containing
     * indices of nearest neighbors and squared distances to said neighbors.
     * @throw Assertion fails if there are not enough points in the tree.
     */
    std::pair<mm::Range<int>, mm::Range<double>> query(const vec_t& point, int k = 1);

    /// Output basic info about given tree.
    template <class V>
    friend std::ostream& operator<<(std::ostream& os, const KDTreeMutable<V>& tree);
};

}  // namespace mm

#endif  // MEDUSA_BITS_SPATIAL_SEARCH_KDTREEMUTABLE_FWD_HPP_
