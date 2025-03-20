#ifndef MEDUSA_BITS_UTILS_PRINT_HPP_
#define MEDUSA_BITS_UTILS_PRINT_HPP_

/**
 * @file
 * Printing helpers for std types.
 */

#include <iostream>
#include <vector>
#include <array>
#include <utility>
#include <tuple>

// additional ostream operators
namespace std {

/// Output pairs as `(1, 2)`. @ingroup utils
template<class T, class U>
std::ostream& operator<<(std::ostream& xx, const std::pair<T, U>& par) {
    return xx << "(" << par.first << "," << par.second << ")";
}

/// Output arrays as `[1, 2, 3]`. @ingroup utils
template<class T, size_t N>
std::ostream& operator<<(std::ostream& xx, const std::array<T, N>& arr) {
    xx << "[";
    for (size_t i = 0; i < N; ++i) {
        xx << arr[i];
        if (i < N - 1) xx << ", ";
    }
    xx << "]";
    return xx;
}

/// Output vectors as `[1, 2, 3]`. @ingroup utils
template<class T, class A>
std::ostream& operator<<(std::ostream& xx, const std::vector<T, A>& arr) {
    // do it like the matlab does it.
    xx << "[";
    for (size_t i = 0; i < arr.size(); ++i) {
        xx << arr[i];
        if (i < arr.size() - 1) xx << ", ";
    }
    xx << "]";
    return xx;
}

/// Output nested vectors as `[[1, 2]; [3, 4]]`. @ingroup utils
template<class T, class A>
std::ostream& operator<<(std::ostream& xx, const std::vector<std::vector<T, A>>& arr) {
    xx << "[";
    for (size_t i = 0; i < arr.size(); ++i) {
        for (size_t j = 0; j < arr[i].size(); ++j) {
            xx << arr[i][j];
            if (j < arr[i].size() - 1) xx << ", ";
        }
        if (i < arr.size() - 1) xx << "; ";
    }
    xx << "]";
    return xx;
}

/// @cond
namespace tuple_print_internal {
template <class Tuple, std::size_t N>
struct TuplePrinter {
    static void print(std::ostream& os, const Tuple& t) {  // recursive
        TuplePrinter<Tuple, N - 1>::print(os, t);
        os << ", " << std::get<N - 1>(t);
    }
};

template <class Tuple>
struct TuplePrinter<Tuple, 1> {
    static void print(std::ostream& os, const Tuple& t) {  // one element
        os << std::get<0>(t);
    }
};

template <class Tuple>
struct TuplePrinter<Tuple, 0> { static void print(std::ostream&, const Tuple&) {} };  // zero elt

}  // namespace tuple_print_internal
/// @endcond

/// Print a tuple as (1, 4.5, abc). @ingroup utils
template <class... Args>
std::ostream& operator<<(std::ostream& os, const std::tuple<Args...>& t) {
    os << "(";
    tuple_print_internal::TuplePrinter<decltype(t), sizeof...(Args)>::print(os, t);
    return os << ")";
}


}  // namespace std

#endif  // MEDUSA_BITS_UTILS_PRINT_HPP_
