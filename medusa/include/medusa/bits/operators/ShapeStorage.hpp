#ifndef MEDUSA_BITS_OPERATORS_SHAPESTORAGE_HPP_
#define MEDUSA_BITS_OPERATORS_SHAPESTORAGE_HPP_

/**
 * @file
 * Implementations of the ShapeStorage interface.
 */

#include "ShapeStorage_fwd.hpp"
#include "printShapes.hpp"

namespace mm {

template <typename Derived, typename vec_t, typename OpFamilies>
int ShapeStorage<Derived, vec_t, OpFamilies>::support(int node, int j) const {
    assert_msg(0 <= node && node <= domain_size_, "Node index %d out or range [0, %d).",
               node, domain_size_);
    assert_msg(0 <= j && j < supportSize(node), "Support index %d out of range [0, %d).",
               j, supportSize(node));
    assert_msg(access(support_, node)[j] != -1, "Attempted access to shapes and support for node "
                                                "%d, which were not computed.", node);
    return access(support_, node)[j];
}

template <typename Derived, typename vec_t, typename OpFamilies>
Eigen::Map<const Eigen::Matrix<int, Eigen::Dynamic, 1>>
ShapeStorage<Derived, vec_t, OpFamilies>::support(int node) const {
    assert_msg(0 <= node && node <= domain_size_, "Node index %d out or range [0, %d).",
               node, domain_size_);
    assert_msg(*access(support_, node) != -1, "Attempted access to shapes and support for node "
                                              "%d, which were not computed.", node);
    return {access(support_, node), supportSize(node)};
}

template <typename Derived, typename vec_t, typename OpFamilies>
Range<int> ShapeStorage<Derived, vec_t, OpFamilies>::supportSizes() const {
    Range<int> sizes(domain_size_, 0);
    for (int i = 0; i < domain_size_; ++i) {
        sizes[i] = supportSize(i);
    }
    return sizes;
}

template <typename Derived, typename vec_t, typename OpFamilies>
Range<int> ShapeStorage<Derived, vec_t, OpFamilies>::supportSizesVec() const {
    Range<int> sizes = supportSizes();
    Range<int> res;
    for (int d = 0; d < dim; ++d) res.append(sizes);
    for (int& i : res) i *= dim;
    return res;
}

template <typename Derived, typename vec_t, typename OpFamilies>
typename ShapeStorage<Derived, vec_t, OpFamilies>::scalar_t
ShapeStorage<Derived, vec_t, OpFamilies>::laplace(int node, int j) const {
    assert_msg(0 <= node && node <= domain_size_, "Node index %d out or range [0, %d).",
               node, domain_size_);
    assert_msg(0 <= j && j < supportSize(node), "Support index %d out of range [0, %d).",
               j, supportSize(node));
    return get<Lap<vec_t::dim>>(node, j);
}

template <typename Derived, typename vec_t, typename OpFamilies>
typename ShapeStorage<Derived, vec_t, OpFamilies>::scalar_t
ShapeStorage<Derived, vec_t, OpFamilies>::d1(int var, int node, int j) const {
    assert_msg(0 <= node && node <= domain_size_, "Node index %d out or range [0, %d).",
               node, domain_size_);
    assert_msg(0 <= j && j < supportSize(node), "Support index %d out of range [0, %d).",
               j, supportSize(node));
    assert_msg(0 <= var && var < dim, "Variable index %d out of bounds [%d, %d).", var, 0, dim);
    return get<Der1s<dim>>(Der1s<dim>::index(Der1<dim>(var)), node, j);
}

template <typename Derived, typename vec_t, typename OpFamilies>
typename ShapeStorage<Derived, vec_t, OpFamilies>::scalar_t
ShapeStorage<Derived, vec_t, OpFamilies>::d2(
        int varmin, int varmax, int node, int j) const {
    assert_msg(0 <= node && node <= domain_size_, "Node index %d out or range [0, %d).",
               node, domain_size_);
    assert_msg(0 <= j && j < supportSize(node), "Support index %d out of range [0, %d).",
               j, supportSize(node));
    assert_msg(0 <= varmin && varmin < dim,
               "Variable varmin %d out of bounds [%d, %d).", varmin, 0, dim);
    assert_msg(0 <= varmax && varmax < dim,
               "Variable varmax %d out of bounds [%d, %d).", varmax, 0, dim);
    assert_msg(varmin <= varmax, "Varmin (%d) must be smaller than varmax (%d).",
               varmin, varmax);
    return get<Der2s<dim>>(Der2s<dim>::index(Der2<dim>(varmin, varmax)), node, j);
}

template <typename Derived, typename vec_t, typename OpFamilies>
Eigen::Map<const Eigen::Matrix<typename ShapeStorage<Derived, vec_t, OpFamilies>::scalar_t,
        Eigen::Dynamic, 1>>
ShapeStorage<Derived, vec_t, OpFamilies>::laplace(int node) const {
    assert_msg(0 <= node && node <= domain_size_, "Node index %d out or range [0, %d).",
               node, domain_size_);
    return getShape<Lap<vec_t::dim>>(node);
}

template <typename Derived, typename vec_t, typename OpFamilies>
Eigen::Map<const Eigen::Matrix<typename ShapeStorage<Derived, vec_t, OpFamilies>::scalar_t,
        Eigen::Dynamic, 1>>
ShapeStorage<Derived, vec_t, OpFamilies>::d1(int var, int node) const {
    assert_msg(0 <= node && node <= domain_size_, "Node index %d out or range [0, %d).",
               node, domain_size_);
    assert_msg(0 <= var && var < dim, "Variable index %d out of bounds [%d, %d).", var, 0, dim);
    return getShape<Der1s<dim>>(Der1s<dim>::index(Der1<dim>(var)), node);
}

template <typename Derived, typename vec_t, typename OpFamilies>
Eigen::Map<const Eigen::Matrix<typename ShapeStorage<Derived, vec_t, OpFamilies>::scalar_t,
        Eigen::Dynamic, 1>>
ShapeStorage<Derived, vec_t, OpFamilies>::d2(int varmin, int varmax, int node) const {
    assert_msg(0 <= node && node <= domain_size_, "Node index %d out or range [0, %d).",
               node, domain_size_);
    assert_msg(0 <= varmin && varmin < dim,
               "Variable varmin %d out of bounds [%d, %d).", varmin, 0, dim);
    assert_msg(0 <= varmax && varmax < dim,
               "Variable varmax %d out of bounds [%d, %d).", varmax, 0, dim);
    assert_msg(varmin <= varmax, "Varmin (%d) must be smaller than varmax (%d).",
               varmin, varmax);
    return getShape<Der2s<dim>>(Der2s<dim>::index(Der2<dim>(varmin, varmax)), node);
}

template <typename Derived, typename vec_t, typename OpFamilies>
void ShapeStorage<Derived, vec_t, OpFamilies>::setSupport(
        int node, const std::vector<int>& support) {
    assert_msg(0 <= node && node <= domain_size_, "Node index %d out or range [0, %d).",
               node, domain_size_);
    assert_msg(static_cast<int>(support.size()) == supportSize(node), "Support of unexpected "
               "length %d, expected %d.", support.size(), supportSize(node));
    std::memcpy(access(support_, node), support.data(),
                support.size()*sizeof(int));
}

template <typename Derived, typename vec_t, typename OpFamilies>
void ShapeStorage<Derived, vec_t, OpFamilies>::setLaplace(
        int node, const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>& shape) {
    assert_msg(0 <= node && node <= domain_size_, "Node index %d out or range [0, %d).",
               node, domain_size_);
    assert_msg(shape.size() == supportSize(node), "Laplace shape of unexpected length %d, "
                                                  "expected %d.", shape.size(), supportSize(node));
    setShape<Lap<vec_t::dim>>(node, shape);
}

template <typename Derived, typename vec_t, typename OpFamilies>
void ShapeStorage<Derived, vec_t, OpFamilies>::setD1(
        int var, int node, const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>& shape) {
    assert_msg(0 <= node && node <= domain_size_, "Node index %d out or range [0, %d).",
               node, domain_size_);
    assert_msg(0 <= var && var < dim, "Variable index %d out of bounds [%d, %d).", var, 0, dim);
    assert_msg(shape.size() == supportSize(node), "D1 shape of unexpected length %d, expected %d.",
               shape.size(), supportSize(node));
    setShape<Der1s<dim>>(var, node, shape);
}

template <typename Derived, typename vec_t, typename OpFamilies>
void ShapeStorage<Derived, vec_t, OpFamilies>::setD2(
        int varmin, int varmax, int node, const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>& shape) {
    assert_msg(0 <= node && node <= domain_size_, "Node index %d out or range [0, %d).",
               node, domain_size_);
    assert_msg(0 <= varmin && varmin < dim,
               "Variable varmin %d out of bounds [%d, %d).", varmin, 0, dim);
    assert_msg(0 <= varmax && varmax < dim,
               "Variable varmax %d out of bounds [%d, %d).", varmax, 0, dim);
    assert_msg(varmin <= varmax, "Varmin (%d) must be smaller than varmax (%d).",
               varmin, varmax);
    assert_msg(shape.size() == supportSize(node), "D2 shape of unexpected length %d, expected %d.",
               shape.size(), supportSize(node));
    int idx = varmax*(varmax+1)/2 + varmin;
    setShape<Der2s<dim>>(idx, node, shape);
}

/// Output basic info about this shape storage.
template <typename D, typename V, typename O>
std::ostream& operator<<(std::ostream& os, const ShapeStorage<D, V, O>& shapes) {
    os << "Shape storage:\n";
    return shapes_internal::printShapes(shapes, os);
}

template <typename Derived, typename vec_t, typename OpFamilies>
size_t ShapeStorage<Derived, vec_t, OpFamilies>::memoryUsed() const {
    size_t size = sizeof(*this);
    for (int i = 0; i < num_operators; ++i) {
        size += mem_used(shapes_[i]);
    }
    return size;
}

/// Contains a helper for resizing storage.
namespace resize_internal {

/// @cond
// A helper for resizing storage -- recursive case
template <typename OpFamilies, size_t i>
struct resize_storage_ {
    template <typename Shapes>
    static void resize(Shapes& shapes, int base_size) {
        resize_storage_<OpFamilies, i-1>::template resize<Shapes>(shapes, base_size);
        std::get<i-1>(shapes).resize(
                base_size*std::tuple_element<i-1, OpFamilies>::type::size(), NaN);
    }
};

// A helper for resizing storage -- base case
template <typename OpFamilies>
struct resize_storage_<OpFamilies, 0> {
    template <typename Shapes> static void resize(Shapes&, int) {}
};
/// @endcond

}  // namespace resize_internal

}  // namespace mm

#endif  // MEDUSA_BITS_OPERATORS_SHAPESTORAGE_HPP_
