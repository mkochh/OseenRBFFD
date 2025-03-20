#ifndef MEDUSA_BITS_OPERATORS_RAGGEDSHAPESTORAGE_HPP_
#define MEDUSA_BITS_OPERATORS_RAGGEDSHAPESTORAGE_HPP_

/**
 * @file
 * Implementation of shape storage for shapes of different sizes.
 */

#include "RaggedShapeStorage_fwd.hpp"
#include "shape_flags.hpp"
#include "ShapeStorage.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <cstring>
#include <Eigen/Core>

namespace mm {

template <typename vec_t, typename OpFamilies>
void RaggedShapeStorage<vec_t, OpFamilies>::resize_(const std::vector<int>& support_sizes) {
    domain_size_ = support_sizes.size();
    support_sizes_ = support_sizes;

    support_starts_.resize(domain_size_, 0);
    for (int i = 1; i < domain_size_; ++i) {
        support_starts_[i] = support_starts_[i-1] + support_sizes_[i-1];
    }
    total_size_ = support_starts_[domain_size_-1] + support_sizes_[domain_size_-1];

    // Fills support with -1, to indicate which nodes have their shapes computed.
    support_.resize(total_size_, -1);
    resize_internal::resize_storage_<OpFamilies, num_operators>::resize(
            shapes_, total_size_);
}

}  // namespace mm

#endif  // MEDUSA_BITS_OPERATORS_RAGGEDSHAPESTORAGE_HPP_
