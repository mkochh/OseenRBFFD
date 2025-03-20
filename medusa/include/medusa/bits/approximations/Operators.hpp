#ifndef MEDUSA_BITS_APPROXIMATIONS_OPERATORS_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_OPERATORS_HPP_

#include "Operators_fwd.hpp"
#include "Monomials_fwd.hpp"
#include "RBFBasis_fwd.hpp"

/**
 * @file
 * Implementations of differential operator families.
 */

namespace mm {
template <typename Derived>
template <typename basis_t>
typename basis_t::scalar_t Operator<Derived>::apply(
        const basis_t&, int, typename basis_t::vector_t,
        const std::vector<typename basis_t::vector_t>&, typename basis_t::scalar_t) const {
    static_assert(!std::is_same<Derived, Derived>::value,
                  "Applying this operator to this basis not implemented.");
}

template <typename Derived>
template <typename basis_t>
typename basis_t::scalar_t Operator<Derived>::applyAt0(
        const basis_t& basis, int index, const std::vector<typename basis_t::vector_t>& support,
        typename basis_t::scalar_t scale) const {
    return apply(basis, typename basis_t::vector_t(0.0), index, support, scale);
}

template <int dimension>
std::array<Der1<dimension>, dimension> Der1s<dimension>::operators() {
    std::array<Der1<dim>, dim> ret;
    for (int d = 0; d < dim; ++d) {
        ret[d] = {d};
    }
    return ret;
}
template <int dimension>
std::array<Der2<dimension>, dimension*(dimension+1)/2> Der2s<dimension>::operators() {
    std::array<Der2<dim>, dim*(dim+1)/2> ret;
    int c = 0;
    for (int d2 = 0; d2 < dim; ++d2) {
        for (int d1 = 0; d1 <= d2; ++d1) {
            ret[c] = {d1, d2};
            ++c;
        }
    }
    return ret;
}

/// Implementation details for analytical operators
namespace operators_internal {
/**
 * Convert variable index to its respective letter, 0 -> x, 1 -> y, etc.
 * For numbers more than 4, generic name x_i is printed.
 */
inline std::string idx_to_letter(int var) {
    if (var < 0 || var >= 4) return "x_" + std::to_string(var);
    return std::string(1, "xyzw"[var]);
}
}  // namespace operators_internal


template <int dimension>
Der1<dimension>::Der1(int var)  : var(var) {
    assert_msg(0 <= var && var <= dim, "Variable %d should be in range [0, %d).", var, dim);
}

template <int dimension>
std::string Der1<dimension>::name() const {
    return format("Der1<%d> wrt. %s", dim, operators_internal::idx_to_letter(var));
}

template <int dimension>
Der2<dimension>::Der2(int var1, int var2) : var1(var1), var2(var2) {
    assert_msg(0 <= var1 && var1 <= dim, "Variable %d should be in range [0, %d).", var1, dim);
    assert_msg(0 <= var2 && var2 <= dim, "Variable %d should be in range [0, %d).", var2, dim);
    assert_msg(var2 >= var1,
            "Second value provided %d should be higher or equal than the first %d", var2, var1);
}

template <int dimension>
std::string Der2<dimension>::name() const {
    return format("Der2<%d> wrt. %s and %s",
            dim, operators_internal::idx_to_letter(var1), operators_internal::idx_to_letter(var2));
}

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_OPERATORS_HPP_
