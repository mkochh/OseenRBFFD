#ifndef MEDUSA_BITS_OPERATORS_SHAPE_FLAGS_HPP_
#define MEDUSA_BITS_OPERATORS_SHAPE_FLAGS_HPP_

/**
 * @file
 * Definition of flags related to shape computation.
 */

#include <string>
#include <tuple>
#include <medusa/bits/approximations/Operators.hpp>

namespace mm {

/// Namespace holding masks for shapes.
namespace sh {

/**
 * Type representing flags for shape functions. It has to be an integral type and combinable
 * with bitwise or operator `|` indicating flag union.
 */
typedef unsigned int shape_flags;
static const shape_flags d1 = 1;  ///< Indicates to calculate d1 shapes.
static const shape_flags lap = 2;  ///< Indicates to calculate laplace shapes.
static const shape_flags d2 = 4;  ///< Indicates to calculate d2 shapes.
static const shape_flags div = d1;  ///< Indicates to prepare all shapes needed for div.
static const shape_flags grad = d1;  ///< Indicates to prepare all shapes needed for grad.
static const shape_flags graddiv = d2;  ///< Indicates to prepare all shapes needed for graddiv.
static const shape_flags all = d1 | d2 | lap;  ///< Indicates to prepare all shapes, default.

/// Convert shape flags to a string representation.
inline std::string str(shape_flags f) {
    std::string s = "";
    if (f & sh::d1) s += " d1";
    if (f & sh::lap) s += " lap";
    if (f & sh::d2) s += " d2";
    return s;
}

/// Converts shape mask to operator tuple type.
template <shape_flags mask, int dim>
struct operator_tuple {
    static_assert(mask < 0, "Type for this mask not defined");
    typedef void type;   ///< The type of the operator tuple.
};
/// @cond
template <int dim> struct operator_tuple<0, dim> {typedef std::tuple<> type; };  // NOLINT
template <int dim> struct operator_tuple<lap, dim> {typedef std::tuple<Lap<dim>> type; };  // NOLINT
template <int dim> struct operator_tuple<lap|d1, dim> { typedef std::tuple<Lap<dim>, Der1s<dim>> type; };  // NOLINT
template <int dim> struct operator_tuple<lap|d2, dim> { typedef std::tuple<Lap<dim>, Der2s<dim>> type; };  // NOLINT
template <int dim> struct operator_tuple<lap|d1|d2, dim> { typedef std::tuple<Lap<dim>, Der1s<dim>, Der2s<dim>> type; };  // NOLINT
template <int dim> struct operator_tuple<d1, dim> { typedef std::tuple<Der1s<dim>> type; };  // NOLINT
template <int dim> struct operator_tuple<d1|d2, dim> { typedef std::tuple<Der1s<dim>, Der2s<dim>> type; };  // NOLINT
template <int dim> struct operator_tuple<d2, dim> { typedef std::tuple<Der2s<dim>> type; };  // NOLINT
/// @endcond

}  // namespace sh

}  // namespace mm

#endif  // MEDUSA_BITS_OPERATORS_SHAPE_FLAGS_HPP_
