#ifndef MEDUSA_BITS_OPERATORS_IMPLICITVECTOROPERATORS_HPP_
#define MEDUSA_BITS_OPERATORS_IMPLICITVECTOROPERATORS_HPP_

/**
 * @file
 * Implementations of implicit vector operators.
 */

#include "ImplicitOperators.hpp"
#include "ImplicitVectorOperators_fwd.hpp"
#include "ShapeStorage_fwd.hpp"

namespace mm {

template <class shape_storage_type, class matrix_type, class rhs_type>
ImplicitVectorOperators<shape_storage_type, matrix_type, rhs_type>::ImplicitVectorOperators(
        const shape_storage_t& ss, matrix_t& M, rhs_t& rhs, int row_offset, int col_offset) :
        ss(&ss), M(&M), rhs(&rhs), row_offset(row_offset), col_offset(col_offset) {
    assert_msg(0 <= row_offset, "Row offset cannot be negative, got %d.", row_offset);
    assert_msg(0 <= col_offset, "Col offset cannot be negative, got %d.", col_offset);
    assert_msg(M.rows() >= dim*ss.size() + row_offset,
               "Matrix does not have enough rows, expected at least %d = %d + %d, got %d.",
               dim*ss.size()+row_offset, row_offset, dim*ss.size(), M.rows());
    assert_msg(M.cols() >= dim*ss.size(),
               "Matrix does not have enough columns, expected at least %d = %d + %d, got %d.",
               dim*ss.size()+col_offset, col_offset, dim*ss.size(), M.cols());
    assert_msg(rhs.size() >= dim*ss.size() + row_offset,
               "RHS vector does not have enough rows, expected at least %d = %d + %d, got %d.",
               dim*ss.size()+row_offset, row_offset, dim*ss.size(), M.rows());
    #ifndef NDEBUG
    // Warning if the matrix or rhs is not zero initialized.
    if (M.squaredNorm() > 1e-10) {
        std::cerr << "Warning: matrix in implicit operators not initialized to zero!" << std::endl;
    }
    if (rhs.squaredNorm() > 1e-10) {
        std::cerr << "Warning: rhs in implicit operators not initialized to zero!" << std::endl;
    }
    #endif
}

/// Output basic info about given operators.
template <typename S, typename M, typename R>
std::ostream& operator<<(std::ostream& os, const ImplicitVectorOperators<S, M, R>& op) {
    os << "Implicit vector operators ";
    os << ((op.hasShapes()) ? "with" : "without") << " linked storage ";
    if (!op.hasShapes()) {
        os << "without linked storage";
    } else {
        os << "over storage: " << *op.ss;
    }
    os << " with " << ((op.hasMatrix()) ? "defined" : "undefined") << " problem matrix";
    os << " and with " << ((op.hasRhs()) ? "defined" : "undefined") << " problem rhs.";
    return os;
}

/// @cond
template <typename Derived, typename vec_t, typename OpFamilies>
template <typename M, typename R>
ImplicitVectorOperators<Derived, M, R>
ShapeStorage<Derived, vec_t, OpFamilies>::implicitVectorOperators(M& matrix, R& rhs) const {
    return {*static_cast<const Derived*>(this), matrix, rhs};
}
/// @endcond
}  // namespace mm

#endif  // MEDUSA_BITS_OPERATORS_IMPLICITVECTOROPERATORS_HPP_
