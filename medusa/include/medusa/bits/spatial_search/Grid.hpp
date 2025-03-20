#ifndef MEDUSA_BITS_SPATIAL_SEARCH_GRID_HPP_
#define MEDUSA_BITS_SPATIAL_SEARCH_GRID_HPP_

/**
 * @file
 * Implementation of Grid.
 */

#include "Grid_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/utils/memutils.hpp>
#include <ostream>

namespace mm {

template <typename T, int dimension, typename IndexType, typename IndexArrayT>
IndexType Grid<T, dimension, IndexType, IndexArrayT>::size(int i) const {
    assert_msg(0 <= i && i < dim, "Dimension index %i must be in range [0, %d)", i, dim);
    return sizes_[i];
}

template <typename T, int dimension, typename IndexType, typename IndexArrayT>
const T& Grid<T, dimension, IndexType, IndexArrayT>::operator()(
        const Grid::IndexArray& index) const {
    assert_msg(inBounds(index, sizes_), "Index %s is out of bounds %s.", index, sizes_);
    return data_[linearIndex(index, sizes_)];
}

template <typename T, int dimension, typename IndexType, typename IndexArrayT>
T& Grid<T, dimension, IndexType, IndexArrayT>::operator()(const Grid::IndexArray& index) {
    assert_msg(inBounds(index, sizes_), "Index %s is out of bounds %s.", index, sizes_);
    return data_[linearIndex(index, sizes_)];
}

template <typename T, int dimension, typename IndexType, typename IndexArrayT>
const T& Grid<T, dimension, IndexType, IndexArrayT>::operator[](const Index& index) const {
    assert_msg(0 <= index && index < size_, "Index %d is out of bounds [0, %d).", index, size_);
    return data_[index];
}

template <typename T, int dimension, typename IndexType, typename IndexArrayT>
T& Grid<T, dimension, IndexType, IndexArrayT>::operator[](const Index& index) {
    assert_msg(0 <= index && index < size_, "Index %d is out of bounds [0, %d).", index, size_);
    return data_[index];
}

template <typename T, int dimension, typename IndexType, typename IndexArrayT>
IndexType Grid<T, dimension, IndexType, IndexArrayT>::computeSize(const Grid::IndexArray& sizes) {
    Index result = 1;
    for (int i = 0; i < dim; ++i) {
        result *= sizes[i];
    }
    return result;
}

template <typename T, int dimension, typename IndexType, typename IndexArrayT>
bool Grid<T, dimension, IndexType, IndexArrayT>::inBounds(const Grid::IndexArray& index,
                                             const Grid::IndexArray& bounds) {
    for (int i = 0; i < dim; ++i) {
        if (index[i] < 0 || index[i] >= bounds[i]) return false;
    }
    return true;
}

template <typename T, int dimension, typename IndexType, typename IndexArrayT>
IndexType Grid<T, dimension, IndexType, IndexArrayT>::linearIndex(
        const Grid::IndexArray& index, const Grid::IndexArray& bounds) {
    assert_msg(inBounds(index, bounds), "Index %s is out of bounds %s.", index, bounds);
    Index result = 0;
    for (int i = 0; i < dim; ++i) {
        result *= bounds[i];
        result += index[i];
    }
    return result;
}

template <typename T, int dimension, typename IndexType, typename IndexArrayT>
IndexArrayT Grid<T, dimension, IndexType, IndexArrayT>::multiIndex(
        Index index, const IndexArray& bounds) {
    IndexArray multi_index;
    for (int i = dim-1; i >= 0; --i) {
        multi_index[i] = index % bounds[i];
        index /= bounds[i];
    }
    return multi_index;
}

/// @cond
template <typename T, int dim, typename I, typename IA>
std::ostream& operator<<(std::ostream& os, const Grid<T, dim, I, IA>& grid) {
    os << dim << "D grid with sizes " << grid.sizes_ << " and " << grid.size_ << " elements, ";
    os << "using " << mem2str(mem_used(grid.data())) << " of memory.";
    return os;
}
/// @endcond

}  // namespace mm

#endif  // MEDUSA_BITS_SPATIAL_SEARCH_GRID_HPP_
