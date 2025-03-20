#ifndef MEDUSA_BITS_UTILS_STDTYPESUTILS_HPP_
#define MEDUSA_BITS_UTILS_STDTYPESUTILS_HPP_

/**
 * @file
 * Declaration of utilities for std types.
 *
 * @example test/utils/stdtypesutils_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/utils/assert.hpp>
#include <string>
#include <vector>
#include <tuple>
#include <type_traits>

namespace mm {

/**
 * Splits string by `delim`, returning a vector of tokens (including empty).
 * @code
 * split("a,b,c,", ",");  // returns {"a", "b", "c", ""}
 * @endcode
 * @ingroup utils
 */
std::vector<std::string> split(const std::string& s, const std::string& delim);
/// Overload for `char`. @sa split
std::vector<std::string> split(const std::string& s, char delim);

/**
 * Joins a vector of strings back together.
 * @param parts Vector of strings.
 * @param joiner String to glue the parts with.
 * @ingroup utils
 */
std::string join(const std::vector<std::string>& parts, const std::string& joiner);
/// Overload for `char`. @sa split
std::string join(const std::vector<std::string>& parts, char joiner);

/// Sorts a container inplace. @ingroup utils
template<typename container_t>
container_t& sort(container_t& v) {
    std::sort(v.begin(), v.end());
    return v;
}

/// Sorts a container inplace according to ordering defined by `pred`. @ingroup utils
template<typename container_t, typename predicate_t>
container_t& sort(container_t& v, const predicate_t& pred) {
    std::sort(v.begin(), v.end(), pred);
    return v;
}

/// Returns a sorted copy of container. @ingroup utils
template<typename container_t>
container_t sorted(container_t v) {
    std::sort(v.begin(), v.end());
    return v;
}

/// Returns a sorted copy of container ordered according to `pred`. @ingroup utils
template<typename container_t, typename predicate_t>
container_t sorted(container_t v, const predicate_t& pred) {
    std::sort(v.begin(), v.end(), pred);
    return v;
}

/**
 * Pads a ragged array with given value.
 * Example:
 * @code
 * pad({{1, 2}, {}, {9, 4, 2}}, -1);  // return {{1, 2, -1}, {-1, -1, -1}, {9, 4, 2}}
 * @endcode
 * @ingroup utils
 */
template <typename container_t, typename T>
container_t pad(container_t container, T value) {
    decltype(container.begin()->size()) maxsize = 0;
    for (const auto& x : container) {
        if (x.size() > maxsize) {
            maxsize = x.size();
        }
    }
    for (auto& x : container) {
        for (auto i = x.size(); i < maxsize; ++i) x.push_back(value);
    }
    return container;
}

/**
 * Returns the first index of type `T` in `Tuple`.
 * @tparam T Type to search for.
 * @tparam Tuple Tuple in which to search.
 * @throws Compile time assertion fails if type `T` is not present.
 */
template <typename T, typename Tuple>
struct tuple_index {
    /// Index of `T` in `tuple`.
    static const std::size_t value = -1;
    static_assert(!std::is_same<Tuple, std::tuple<>>::value,  // Did you access invalid operators?
                  "Could not get index of type `T` in given `Tuple`");
};

/// Succesful match.
template <typename T, typename... Types>
struct tuple_index<T, std::tuple<T, Types...>> {
    static const std::size_t value = 0;  ///< Index of `T` in `tuple`.
};

/// Unsuccessful match and recursive search.
template <typename T, typename U, typename... Types>
struct tuple_index<T, std::tuple<U, Types...>> {
    /// Index of `T` in `tuple`.
    static const std::size_t value = 1 + tuple_index<T, std::tuple<Types...>>::value;
};

/// Find type T in `Tuple` -- declaration
template <typename T, typename Tuple>
struct tuple_has_type;

/// Find type T in `Tuple` -- failed case
template <typename T>
struct tuple_has_type<T, std::tuple<>> : std::false_type {};

/// Find type T in `Tuple` -- unsuccessful, continue searching
template <typename T, typename U, typename... Ts>
struct tuple_has_type<T, std::tuple<U, Ts...>> : tuple_has_type<T, std::tuple<Ts...>> {};

/// Find type T in `Tuple` -- success
template <typename T, typename... Ts>
struct tuple_has_type<T, std::tuple<T, Ts...>> : std::true_type {};

}  // namespace mm

#endif  // MEDUSA_BITS_UTILS_STDTYPESUTILS_HPP_
