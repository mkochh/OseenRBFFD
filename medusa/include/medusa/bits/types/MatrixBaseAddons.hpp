#ifndef MEDUSA_BITS_TYPES_MATRIXBASEADDONS_HPP_
#define MEDUSA_BITS_TYPES_MATRIXBASEADDONS_HPP_

/**
 * @file
 * Plugins for `Eigen::MatrixBase` class.
 *
 * This extensions add lexicographical compare to vectors, so that they can be sorted.
 * They also add comparisons to scalars, returning the set of indices, for which the comparison
 * holds. Some member typedefs are added as well.
 */

enum { /** Number of elements of this matrix. Negative if dynamic */
    dim = RowsAtCompileTime * ColsAtCompileTime };

typedef Scalar scalar_t;   ///< Type of the elements, alias of `Scalar`.

/// Lexicographical compare of vectors.
template<typename OtherDerived>
bool operator<(const MatrixBase<OtherDerived>& arg) const {
    for (int i = 0; i < this->size(); ++i) {
        if (this->operator[](i) == arg[i]) continue;
        return this->operator[](i) < arg[i];
    }
    return false;
}

/// Lexicographical compare of vectors.
template<typename OtherDerived>
bool operator>(const MatrixBase<OtherDerived>& arg) const { return arg < *this; }

/// Lexicographical compare of vectors.
template<typename OtherDerived>
bool operator<=(const MatrixBase<OtherDerived>& arg) const { return !(arg < *this); }

/// Lexicographical compare of vectors.
template<typename OtherDerived>
bool operator>=(const MatrixBase<OtherDerived>& arg) const { return !(*this < arg); }

/**
 * Returns list of indexes for which predicate returns `true`.
 * @tparam Pred Any callable object, e.g.\ lambda func, functional,
 * class with `operator()` defined.
 */
template<class Pred>
mm::indexes_t filter(const Pred& pred) const {
    mm::indexes_t ret;
    for (int i = 0; i < size(); ++i)
        if (pred(this->operator[](i))) ret.push_back(i);
    return ret;
}

/// Returns list of indexes of elements that are lower than `v`.
mm::indexes_t operator<(const Scalar& v) const {
    mm::indexes_t ret;
    for (int i = 0; i < size(); ++i)
        if (this->operator[](i) < v) ret.push_back(i);
    return ret;
}

/// Returns list of indexes of elements that are greater than `v`.
mm::indexes_t operator>(const Scalar& v) const {
    mm::indexes_t ret;
    for (int i = 0; i < size(); ++i)
        if (this->operator[](i) > v) ret.push_back(i);
    return ret;
}

/// Returns list of indexes of elements that are lower or equal than `v`.
mm::indexes_t operator<=(const Scalar& v) const {
    mm::indexes_t ret;
    for (int i = 0; i < size(); ++i)
        if (this->operator[](i) <= v) ret.push_back(i);
    return ret;
}

/// Returns list of indexes of elements that are greater or equal than `v`.
mm::indexes_t operator>=(const Scalar& v) const {
    mm::indexes_t ret;
    for (int i = 0; i < size(); ++i)
        if (this->operator[](i) >= v) ret.push_back(i);
    return ret;
}

/// Returns list of indexes for which their elements compare equal to a
mm::indexes_t operator==(const Scalar& v) const {
    mm::indexes_t ret;
    for (int i = 0; i < size(); ++i)
        if (this->operator[](i) == v) ret.push_back(i);
    return ret;
}

/// Returns list of indexes for which their elements compare not equal to a
mm::indexes_t operator!=(const Scalar& v) const {
    mm::indexes_t ret;
    for (int i = 0; i < size(); ++i)
        if (this->operator[](i) != v) ret.push_back(i);
    return ret;
}

#endif  // MEDUSA_BITS_TYPES_MATRIXBASEADDONS_HPP_
