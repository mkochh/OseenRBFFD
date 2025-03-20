#ifndef MEDUSA_BITS_TYPES_VECTORFIELD_FWD_HPP_
#define MEDUSA_BITS_TYPES_VECTORFIELD_FWD_HPP_

/**
 * @file
 * Declaration of VectorField.
 *
 * @example test/types/VectorField_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/utils/assert.hpp>
#include "Vec_fwd.hpp"
#include "traits.hpp"
#include <Eigen/Core>

namespace mm {

/**
 * Represents a discretization of a vector field, a finite collection of vectors.
 * The vector field is stored in a ColMajor fashion.
 * This type is fully compatible with Eigen types and can use the same API.
 * It also supports iteration over columns.
 *
 * @tparam Scalar
 * @tparam dimension Dimensionality of the vector field.
 *
 * Usage example:
 * @snippet types/VectorField_test.cpp Vector field usage example
 * @ingroup types
 *
 * @warning
 * The `.size()` method inherited from Eigen returns the number of all elements,
 * e. g. number of vectors multiplied by the dimension.
 * Use `.rows()` method to get only the number of vectors in the vector field.
 */
template <typename Scalar, int dimension>
class VectorField : public Eigen::Matrix<Scalar, Eigen::Dynamic, dimension> {
  public:
    typedef Scalar scalar_t;  ///< Scalar data type.
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, dimension> Base;  ///< Base class.
    using Base::operator=;  // inherit assignment operators

    /// Construct a vector field of size 0.
    VectorField(void) : Base() {}
    /// Default copy constructor.
    VectorField(const VectorField&) = default;
    /// Default copy assignment.
    VectorField& operator=(const VectorField&) = default;
    /// Default move constructor.
    VectorField(VectorField&&) noexcept = default;
    /// Default move assignment.
    VectorField& operator=(VectorField&&) noexcept = default;

    /// Allows construction of VectorField from Eigen expressions.
    template <typename OtherDerived>
    VectorField(const Eigen::MatrixBase<OtherDerived>& other);

    /// Assignments of Eigen expressions to VectorField.
    template <typename OtherDerived>
    VectorField& operator=(const Eigen::MatrixBase<OtherDerived>& other);

    /**
     * Construct a vector field on `N` points. Allocates the appropriate space,
     * the values are initialized as specified by Eigen's default initialization.
     */
    template <typename SizeType>
    explicit VectorField(SizeType N) : Base(N, dimension) {}

    /// Assigns vector field to be a constant vector `v`.
    VectorField& operator=(const Vec<Scalar, dimension>& v);

    /// Construct VectorField from its linear representation, obtained from asLinear().
    template <typename OtherDerived>
    static VectorField fromLinear(const Eigen::PlainObjectBase<OtherDerived>& other);

  private:
    /// Represents a non contiguous view to a vector field.
    class VectorFieldView {
        VectorField& sf;  ///< Reference to the viewed field.
        const indexes_t& indexes;  ///< Indexes of this view.
        friend class VectorField;
        /// Construct a non contiguous view to a scalar field.
        VectorFieldView(VectorField& sf, const indexes_t& indexes) : sf(sf), indexes(indexes) {}
      public:
        /// Assign a vector to this view.
        void operator=(const Vec<scalar_t, dimension>& other) {
            sf(indexes, Eigen::all) = other.transpose().replicate(indexes.size(), 1);
        }
    };

  public:
    /// Return a single row as a vector, non-const version.
    Eigen::Transpose<typename Base::RowXpr> operator[](typename Eigen::Index i) {
        return Base::row(i).transpose();
    }
    /// Return a single row as a vector, const version.
    SINL Eigen::Transpose<typename Base::ConstRowXpr> operator[](typename Eigen::Index i) const {
        return Base::row(i).transpose();
    }

    using Base::operator();
    /// Return a single row as a vector, non-const version.
    SINL Eigen::Transpose<typename Base::RowXpr> operator()(typename Eigen::Index i) {
        return Base::row(i).transpose();
    }
    /// Return a single row as a vector, const version.
    SINL Eigen::Transpose<typename Base::ConstRowXpr> operator()(typename Eigen::Index i) const {
        return Base::row(i).transpose();
    }

    /**
     * Return the `i`-th component of this vector field. The result of this function can be
     * treated as a scalar field.
     */
    SINL typename Base::ColXpr c(typename Eigen::Index i) { return Base::col(i); }
    /// Const version of VectorField::c(). @sa c
    SINL typename Base::ConstColXpr c(typename Eigen::Index i) const { return Base::col(i); }

    /// Return an indexed view to a subset of vectors.
    Eigen::IndexedView<Base, indexes_t, Eigen::internal::AllRange<dimension>>
    SINL operator()(const indexes_t& rowIndices) {
        return Base::operator()(rowIndices, Eigen::all);
    }

    /// Returns an indexed view to a subset of vectors, that can be assigned.
    SINL VectorFieldView operator[](const indexes_t& indexes) {
        return VectorFieldView(*this, indexes);
    }

    /// Pointer to the first coefficient.
    Scalar* begin() { return this->data(); }
    /// Pointer to the first coefficient.
    const Scalar* begin() const { return this->data(); }

    /// Pointer to the last coefficient.
    Scalar* end() { return this->data() + this->size(); }
    /// Const pointer to the last coefficient.
    const Scalar* end() const { return this->data() + this->size(); }

    /// Const version of @ref begin.
    const Scalar* cbegin() const { return this->data(); }
    /// Const version of @ref end.
    const Scalar* cend() const { return this->data() + this->size(); }
};
/// The vector_type trait definition for VectorField.
template <typename Scalar, int dim> struct vector_type<VectorField<Scalar, dim>> {
    typedef Vec<Scalar, dim> type;  ///< Underlying vector type.
};

typedef VectorField<double, 1> VectorField1d;  ///< One dimensional vector field of doubles.
typedef VectorField<double, 2> VectorField2d;  ///< Two dimensional vector field of doubles.
typedef VectorField<double, 3> VectorField3d;  ///< Three dimensional vector field of doubles.

}  // namespace mm

#endif  // MEDUSA_BITS_TYPES_VECTORFIELD_FWD_HPP_
