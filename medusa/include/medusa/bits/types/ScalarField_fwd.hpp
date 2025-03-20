#ifndef MEDUSA_BITS_TYPES_SCALARFIELD_FWD_HPP_
#define MEDUSA_BITS_TYPES_SCALARFIELD_FWD_HPP_

/**
 * @file
 * Declaration of ScalarField.
 *
 * @example test/types/ScalarField_test.cpp
 */

#include <medusa/Config.hpp>
#include <Eigen/Core>
#include "traits.hpp"

namespace mm {

/**
 * Represents a discretization of a scalar field, a finite collection of scalars.
 * This type is fully compatible with Eigen types and can use the same API.
 * @tparam Scalar Type of elements in the field.
 *
 * Usage example:
 * @snippet types/ScalarField_test.cpp Scalar field usage example
 * @ingroup types
 */
template <typename Scalar>
class ScalarField : public Eigen::Matrix<Scalar, Eigen::Dynamic, 1> {
  public:
    typedef Scalar scalar_t;  ///< Scalar data type.
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Base;  ///< Base class typedef
    using Base::Base;  // inherit constructors
    using Base::operator=;  // inherit assignment operators
    /// This constructor allows construction of ScalarField from Eigen expressions.
    template <typename OtherDerived>
    ScalarField(const Eigen::MatrixBase<OtherDerived>& other);
    /// This method allows assignments of Eigen expressions to ScalarField.
    template <typename OtherDerived>
    ScalarField& operator=(const Eigen::MatrixBase<OtherDerived>& other);

  private:
    /// Represents a non contiguous view to a scalar field.
    class ScalarFieldView {
        ScalarField& sf;  ///< Reference to the viewed field.
        const indexes_t& indexes;  ///< Indexes of this view.
        friend class ScalarField;
        /// Construct a non contiguous view to a scalar field.
        ScalarFieldView(ScalarField& sf, const indexes_t& indexes) : sf(sf), indexes(indexes) {}
      public:
        /// Assign a scalar to this view.
        void operator=(const scalar_t& s) {
            sf(indexes) = ScalarField::Constant(indexes.size(), s);
        }
        /// Assign another Eigen expression to this view.
        template <typename OtherDerived>
        void operator=(const Eigen::MatrixBase<OtherDerived>& other) {
            sf(indexes) = other;
        }
        /// Assign another ScalarFieldView. @warning Aliasing may occur.
        void operator=(const ScalarFieldView& other) {
            sf(indexes) = other.sf(other.indexes);
        }
    };

  public:
    using Base::operator[];
    /// Multiindex view to given ScalarField.
    ScalarFieldView operator[](const indexes_t& indexes) {
        return ScalarFieldView(*this, indexes);
    }
};

typedef ScalarField<double> ScalarFieldd;  ///< Convenience typedef for ScalarField of doubles.
typedef ScalarField<float> ScalarFieldf;  ///< Convenience typedef for ScalarField of floats.

/// The scalar_type trait definition for ScalarField.
template <typename T> struct scalar_type<ScalarField<T>> {
    typedef T type;  ///< Underlying scalar type.
};

}  // namespace mm

#endif  // MEDUSA_BITS_TYPES_SCALARFIELD_FWD_HPP_
