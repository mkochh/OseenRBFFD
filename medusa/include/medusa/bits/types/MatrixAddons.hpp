#ifndef MEDUSA_BITS_TYPES_MATRIXADDONS_HPP_
#define MEDUSA_BITS_TYPES_MATRIXADDONS_HPP_

/**
 * @file
 * Plugins for `Eigen::Matrix` class.
 *
 * This extension adds some convenience constructors and assignments from scalars and
 * initializer lists for fixed size matrices. Additionally if makes all matrices compatible with
 * range based `for` loops and STL containers by adding appropriate `begin()` and `end()` methods.
 * It also offers linear and row linear views to matrices.
 */

/// Assign scalar to a matrix, setting all its elements to `s`.
Matrix& operator=(const Scalar& s) {
    this->setConstant(s);
    return *this;
}

/// Construct matrix from scalar. Enabled only for fixed size matrices.
Matrix(const Scalar& s) {
    _init1alt<RowsAtCompileTime != Dynamic && ColsAtCompileTime != Dynamic>(s);
}

private:
/// @cond
template <bool is_fixed>
void _init1alt(const Scalar& s,
               typename std::enable_if<is_fixed>::type* = nullptr) {
    this->setConstant(s);
}

template <bool is_fixed>
void _init1alt(const Scalar& s, typename std::enable_if<!is_fixed>::type* = nullptr) {
    Base::_check_template_params();
    Base::template _init1<Scalar>(s);
}
/// @endcond

public:
/// Construct a fixed sized matrix from an initializer list.
Matrix(std::initializer_list<Scalar> lst) {
    EIGEN_STATIC_ASSERT_FIXED_SIZE(Matrix);
    assert((lst.size() == RowsAtCompileTime || lst.size() == ColsAtCompileTime) &&
           "Initializer list of inappropriate size.");
    int i = 0;
    for (const Scalar& s : lst) {
        this->operator()(i++) = s;
    }
}

/// Returns a view to a matrix as a column vector.
Map<Matrix<Scalar, Dynamic, 1>> asLinear() {
    return {this->data(), this->size()};
}
Map<Matrix<Scalar, Dynamic, 1>> asLinear() const {
    return {this->data(), this->size()};
}
/// Returns a view to a matrix as a row vector.
Map<Matrix<Scalar, Dynamic, 1>> asRowLinear() {
    return {this->data(), this->size()};
}
Map<Matrix<Scalar, Dynamic, 1>> asRowLinear() const {
    return {this->data(), this->size()};
}

#endif  // MEDUSA_BITS_TYPES_MATRIXADDONS_HPP_
