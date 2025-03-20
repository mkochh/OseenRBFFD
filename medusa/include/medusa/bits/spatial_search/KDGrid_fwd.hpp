#ifndef MEDUSA_BITS_SPATIAL_SEARCH_KDGRID_FWD_HPP_
#define MEDUSA_BITS_SPATIAL_SEARCH_KDGRID_FWD_HPP_

/**
 * @file
 * Declaration of KDGrid.
 *
 * @example test/spatial_search/KDGrid_test.cpp
 */

#include "Grid_fwd.hpp"
#include <iosfwd>
#include <vector>

namespace mm {

/**
 * Search structure over given d-dimensional box with given cell size.
 * At most one point per cell can be stored. This is often a faster and more memory hungry
 * substitute for KDTree.
 *
 * Usage example:
 * @snippet spatial_search/KDGrid_test.cpp KDGrid usage example
 * @ingroup utils
 *
 * @sa GeneralFill, KDTreeMutable
 */
template <typename vec_t>
class KDGrid {
  public:
    typedef vec_t vector_t;   ///< Vector type used.
    typedef typename vector_t::scalar_t scalar_t;   ///< Scalar type used.
    /// Store dimension of the space.
    enum { /** Dimensionality of the space. */ dim = vec_t::dim };

  private:
    vec_t bot_;  ///< Lower bound of the box.
    vec_t top_;  ///< Upper bound of the box.
    scalar_t cell_size_;  ///< Size of a single cell.
    Grid<int, dim> grid_;  ///< Background grid of cells pointing to sequential point indices.
    std::vector<vec_t> points_;  ///< List of inserted points.
    typedef typename Grid<int, dim>::IndexArray IndexArray;   ///< Multi-index type

  public:
    /// Construct a search grid from `bot` to `top` with given `cell_size`.
    KDGrid(const vec_t& bot, const vec_t& top, scalar_t cell_size) :
            bot_(bot), top_(top), cell_size_(cell_size),
            grid_(compute_size(bot, top, cell_size), -1),
            points_() {}

    /// Get bottom bound.
    const vec_t& bottom() const { return bot_; }
    /// Get top bound.
    const vec_t& top() const { return top_; }
    /// Get cell size.
    scalar_t cellSize() const { return cell_size_; }
    /// Get access to underlying grid of point indices.
    const Grid<int, dim>& grid() const { return grid_; }
    /// Get access to array of stored points.
    const std::vector<vec_t>& points() const { return points_; }
    /// Get point by its sequential index.
    const vec_t& point(int idx) const { return points_[idx]; }

    /// Get number of points stored in the structure.
    int size() const { return points_.size(); }

    /**
     * Insert a new point in the structure.
     * @param p The point to insert.
     * @return Sequential index of the point.
     * @throw Assertion fails if there is also point in this cell or if the point is out of range.
     */
    int insert(const vec_t& p) { return insert(p, compute_index(p)); }

    /// Vectorised version of @ref insert
    void insert(const std::vector<vec_t>& pts) {
        for (const vec_t& p : pts) insert(p);
    }

    /// Check if any point exists in sphere centered at `p` with radius `r`.
    bool existsPointInSphere(const vec_t& p, scalar_t r) {
        return existsPointInSphere(p, r, compute_index(p));
    }

    /// Output some information about the search grid.
    template <typename V>
    friend std::ostream& operator<<(std::ostream& os, const KDGrid<V>& search);

  private:
    /// Compute the sizes of the rid along each dimension.
    static IndexArray compute_size(const vec_t& bot, const vec_t& top, scalar_t cell_size);

    /// Compute point multi-index.
    IndexArray compute_index(const vec_t& p);

    /// Implementation of @ref existsPointInSphere
    bool existsPointInSphere(const vec_t& p, scalar_t r, const IndexArray& index);

    /// Implementation of @ref insert
    int insert(const vec_t& p, const IndexArray& index);

    /// Increment a dim-ary counter.
    bool increment(IndexArray& cur, int span, const IndexArray& ref);
};

}  // namespace mm

#endif  // MEDUSA_BITS_SPATIAL_SEARCH_KDGRID_FWD_HPP_
