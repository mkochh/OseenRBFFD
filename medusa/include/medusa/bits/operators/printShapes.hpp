#ifndef MEDUSA_BITS_OPERATORS_PRINTSHAPES_HPP_
#define MEDUSA_BITS_OPERATORS_PRINTSHAPES_HPP_

/**
 * @file
 * Suppport functions for shape printing.
 */

#include "shape_flags.hpp"
#include <iostream>
#include <medusa/bits/utils/memutils.hpp>

namespace mm {
namespace shapes_internal {

/// @cond
// recursive case
template <typename OpFamilies, size_t i>
struct name_printer {
    static void print(std::ostream& os) {
        name_printer<OpFamilies, i-1>::print(os);
        os << ", " << std::tuple_element<i-1, OpFamilies>::type::name();
    }
};

// base case 1 element
template <typename OpFamilies>
struct name_printer<OpFamilies, 1> {
    static void print(std::ostream& os) {
        os << std::tuple_element<0, OpFamilies>::type::name();
    }
};

// base case 0 element
template <typename OpFamilies>
struct name_printer<OpFamilies, 0> {
    static void print(std::ostream& os) {
        os << "()";
    }
};

// recursive case
template <typename shape_storage_t, size_t op>
struct shape_printer {
    static void print(std::ostream& os, const shape_storage_t& storage) {
        shape_printer<shape_storage_t, op-1>::print(os, storage);
        int N = std::min(storage.size(), 5);
        using el = typename std::tuple_element<
                op-1, typename shape_storage_t::op_families_tuple>::type;
        os << "    " << el::type_name() << " (family #" << op-1 << ") shape sample:\n";
        for (int c = 0; c < el::size(); ++c) {
            os << "      Operator #" << c << ":\n";
            for (int i = 0; i < N; ++i) {
                os << "        node " << i << ":";
                for (int j = 0; j < storage.supportSize(i); ++j) {
                    os << " " << storage.template get<op-1>(c, i, j);
                }
                os << '\n';
            }
        }
    }
};

// base case
template <typename shape_storage_t> struct shape_printer<shape_storage_t, 0> {
    static void print(std::ostream&, const shape_storage_t&) {}
};
/// @endcond

/// Output basic info about given shape storage to `os`.
template <typename shape_storage_t>
std::ostream& printShapes(const shape_storage_t& storage, std::ostream& os) {
    os << "  dimension: " << storage.dim << '\n'
       << "  domain size: " << storage.size() << '\n'
       << "  operator families = ";
    name_printer<typename shape_storage_t::op_families_tuple,
                 shape_storage_t::num_operators>::print(os);
    os << "\n";
    shape_printer<shape_storage_t, shape_storage_t::num_operators>::print(os, storage);

    std::size_t memory_used = storage.memoryUsed();
    os << "  mem used total: " << mem2str(memory_used);
    return os;
}

}  // namespace shapes_internal
}  // namespace mm

#endif  // MEDUSA_BITS_OPERATORS_PRINTSHAPES_HPP_
