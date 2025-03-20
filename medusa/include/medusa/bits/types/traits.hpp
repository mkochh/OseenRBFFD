#ifndef MEDUSA_BITS_TYPES_TRAITS_HPP_
#define MEDUSA_BITS_TYPES_TRAITS_HPP_

/**
 * @file
 * Type traits definition for vector and scalar fields.
 */

namespace mm {

/**
 * Type trait for scalar fields to obtain their underlying scalar type. This trait should be defined
 * for each type you want to use as a scalar field with explicit operators. It is
 * already defined for Range, std::vector and ScalarField.
 */
template <typename scalar_field_t>
struct scalar_type {
    /// Default scalar type.
    typedef typename scalar_field_t::Scalar type;  // Reason for error: Underlying scalar type was not defined for given scalar_field_t. Please specify the required 'scalar_type' type trait.  // NOLINT(*)
};

/**
 * Type trait for vector fields to obtain their underlying vector type. This trait should be defined
 * for each type you want to use as a vector field with explicit operators. It is
 * already defined for Range, std::vector and VectorField.
 */
template <typename vector_field_t>
struct vector_type {
    static_assert(static_cast<vector_field_t*>(nullptr),
                  "Underlying scalar type was not defined for given vector_field_t. Please specify "
                  "the required 'vector_type' type trait.");
    /// By default invalid underlying vector type. Must be specified for each type separately.
    typedef void type;
};

}  // namespace mm

#endif  // MEDUSA_BITS_TYPES_TRAITS_HPP_
