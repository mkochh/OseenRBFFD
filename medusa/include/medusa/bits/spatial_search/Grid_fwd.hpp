#ifndef MEDUSA_BITS_SPATIAL_SEARCH_GRID_FWD_HPP_
#define MEDUSA_BITS_SPATIAL_SEARCH_GRID_FWD_HPP_


/**
 * @file
 * Declaration of Grid.
 *
 * @example test/spatial_search/Grid_test.cpp
 */


#include <medusa/Config.hpp>
#include <array>
#include <vector>
#include <iosfwd>

namespace mm {

/**
 * Class representing a simple n-dimensional grid structure, which supports indexing and storing
 * values. The values are stored in row-major order.
 *
 * @tparam T Type of the values stored in the grid.
 * @tparam dimension Dimension of the grid.
 * @tparam IndexType Type used for indexing.
 *
 * Usage example:
 * @snippet spatial_search/Grid_test.cpp Grid usage example
 * @ingroup utils
 *
 * @sa KDGrid
 */
template <typename T, int dimension, typename IndexType = int,
          typename IndexArrayT = std::array<IndexType, dimension>>
class Grid {
  public:
    typedef T value_type;  ///< Type of the values stored in the grid.
    typedef IndexType Index;  ///< Type of the index used.
    enum { /** Dimensionality of the domain. */ dim = dimension };
    typedef IndexArrayT IndexArray;   ///< Multiindex type.

  private:
    IndexArray sizes_;  ///< Boundary of the grid in all dimensions.
    IndexType size_;  ///< Total number of grid cells.
    std::vector<T> data_;  ///< Data stored in the grid.

  public:
    /// Construct a zero initialized grid with given sizes.
    Grid(const IndexArray& sizes) : sizes_(sizes), size_(computeSize(sizes)), data_(size_) {}
    /// Construct a grid with given sizes initialized to `value`.
    Grid(const IndexArray& sizes, const T& value) :
            sizes_(sizes), size_(computeSize(sizes)), data_(size_, value) {}

    /// Get grid sizes.
    IndexArray sizes() const { return sizes_; }
    /// Get side in i-th dimension.
    inline Index size(int i) const;
    /// Get total number of elements.
    Index size() const { return size_; }

    /// Get read-only access to linear data.
    const std::vector<T>& data() const { return data_; }
    /// Get read-write access to linear data.
    std::vector<T>& data() { return data_; }

    /// Readonly access to the grid. @throw Assertion fails if index if out of bounds.
    inline const T& operator()(const IndexArray& index) const;

    /// Read-write access to the grid. @throw Assertion fails if index if out of bounds.
    inline T& operator()(const IndexArray& index);

    /// Readonly access to the grid. @throw Assertion fails if index if out of bounds.
    inline const T& operator[](const Index& index) const;

    /// Read-write access to the grid. @throw Assertion fails if index if out of bounds.
    inline T& operator[](const Index& index);

    /// Compute linear index from given multi-index. @sa multiIndex
    Index linearIndex(const IndexArray& index) const { return linearIndex(index, sizes_); }
    /// Compute multi-index from given linear index. @sa linearIndex
    IndexArray multiIndex(const Index& index) const { return multiIndex(index, sizes_); }

    /// Check if given index in in bounds.
    bool inBounds(const IndexArray& index) const { return inBounds(index, sizes_); }

  private:
    /// Compute the number of elements.
    static Index computeSize(const IndexArray& sizes);

    /// Check if index is in given bounds.
    static bool inBounds(const IndexArray& index, const IndexArray& bounds);

    /// Compute the linear index with respect to given bounds.
    static Index linearIndex(const IndexArray& index, const IndexArray& bounds);
    /// Compute the multi-index with respect to given bounds.
    static IndexArray multiIndex(Index index, const IndexArray& bounds);

  public:
    /// Output some information about given grid.
    template <typename U, int D, typename I, typename IA>
    friend std::ostream& operator<<(std::ostream& os, const Grid<U, D, I, IA>& grid);
};

}  // namespace mm

#endif  // MEDUSA_BITS_SPATIAL_SEARCH_GRID_FWD_HPP_
