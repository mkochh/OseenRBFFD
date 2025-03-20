#ifndef MEDUSA_BITS_OPERATORS_UNIFORMSHAPESTORAGE_HPP_
#define MEDUSA_BITS_OPERATORS_UNIFORMSHAPESTORAGE_HPP_

/**
 * @file
 * Implementation of shape storage for uniformly long shapes.
 */

#include "UniformShapeStorage_fwd.hpp"
#include "ShapeStorage.hpp"
#include "printShapes.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <cstring>
#include <Eigen/Core>

namespace mm {

template <typename vec_t, typename OpFamilies>
void UniformShapeStorage<vec_t, OpFamilies>::resize_(const std::vector<int>& support_sizes) {
    support_size_ = (domain_size_ > 0) ? support_sizes[0] : 0;
    for (int s : support_sizes) {
        assert_msg(support_size_ == s, "Not all support sizes are equal, got sizes %d and %d. "
                                       "Use RaggedShapeStorage instead.", support_size_, s);
    }

    // Fills support with -1, to indicate which nodes have their shapes computed.
    support_.resize(domain_size_ * support_size_, -1);
    resize_internal::resize_storage_<OpFamilies, num_operators>::resize(
            shapes_, domain_size_ * support_size_);
}

}  // namespace mm

#endif  // MEDUSA_BITS_OPERATORS_UNIFORMSHAPESTORAGE_HPP_
