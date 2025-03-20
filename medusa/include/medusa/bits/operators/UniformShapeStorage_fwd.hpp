#ifndef MEDUSA_BITS_OPERATORS_UNIFORMSHAPESTORAGE_FWD_HPP_
#define MEDUSA_BITS_OPERATORS_UNIFORMSHAPESTORAGE_FWD_HPP_

/**
 * @file
 * Declaration of shape storage for uniformly long shapes.
 *
 * @example test/operators/UniformShapeStorage_test.cpp
 */

#include <medusa/Config.hpp>
#include <vector>
#include <medusa/bits/operators/shape_flags.hpp>
#include <iosfwd>
#include <Eigen/Core>
#include "ShapeStorage_fwd.hpp"

namespace mm {

/**
 * Efficiently stores shape functions of uniform length.
 * This class is used to store shape functions (stencil weights)
 * generated for discretizations where supports of all nodes have equal size
 * e.g.\ the support consists of 9 closest nodes. This class
 * is more efficient than storing the shapes in a nested type, such as
 * `std::vector<std::vector<T>>`, see
 * [technical report](http://e6.ijs.si/medusa/wiki/index.php/File:Tech_report.pdf).
 *
 * If supports of node vary in size, use RaggedShapeStorage instead.
 * *
 * @tparam vec_t Vector type used in computations, specifies the dimensionality of the domain
 * and scalar type for numerical computations.
 * @tparam OpFamilies A list of operator families for which the shapes will be stored.
 * The basic operator families are @ref Lap, @ref Der1s and @ref Der2s representing
 * the Laplacian, 1st and 2nd derivatives, respectively.
 * All these operators are computed if the template parameters are not explicitly specified.
 *
 * If you try to call a function that need other shapes than the computed ones,
 * you will get a *compile time* error like:
 * @code
 * error: static_assert failed due to requirement '!std::is_same<std::tuple<>, std::tuple<> >::value' "Could not find type `T` in given `Tuple`"
 *  static_assert(!std::is_same<Tuple, std::tuple<>>::value,  // Did you access invalid operators?
 *  ^             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * @endcode
 *
 * Usage example:
 * @snippet operators/UniformShapeStorage_test.cpp Uniform shape storage usage example
 * @ingroup operators
 *
 * @sa sh, RaggedShapeStorage, ShapeStorage, computeShapes
 */

template <typename vec_t, typename OpFamilies =
        std::tuple<Lap<vec_t::dim>, Der1s<vec_t::dim>, Der2s<vec_t::dim>>>
class UniformShapeStorage :
        public ShapeStorage<UniformShapeStorage<vec_t, OpFamilies>, vec_t, OpFamilies> {
  public:
    typedef vec_t vector_t;  ///< Vector type used.
    typedef typename vec_t::scalar_t scalar_t;  ///< Scalar type used.
    /// Bitmask telling us which shapes to create.
    enum { /** Dimensionality of the domain. */ dim = vec_t::dim };

  private:
    /// Parent class
    typedef ShapeStorage<UniformShapeStorage<vec_t, OpFamilies>, vec_t, OpFamilies> base_t;
    friend base_t;  ///< Be friends with derived class.
    using base_t::domain_size_;
    using base_t::support_;
    using base_t::shapes_;

    /// Support size.
    int support_size_;

  public:
    using base_t::num_operators;  ///< Number of operators stored in this storage.
    using base_t::size;  ///< Number of nodes that shapes can be stored for.

    /// Constructs an empty shape storage with @ref size 0.
    UniformShapeStorage() : support_size_(0) {}

    /// Returns support size of `node`-th node.
    int supportSize(int /* node */) const { return support_size_; }

  private:
    /**
     * Resizes the storage to accommodate shapes of given sizes. If support sizes are `{9, 9, 9}`
     * the class will allocate space for shapes for 3 nodes with 9 support nodes each.
     * The containers are zero initialized.
     * @throws Assertion fails if the elements of `support_sizes` are not all the same.
     */
    void resize_(const std::vector<int>& support_sizes);

    /// Returns pointer to the start of values for `node`-th node for `op`-th operator.
    template <typename T> T* access(std::vector<T>& v, int op, int node) const {
        return v.data() + op*domain_size_*support_size_ + node*support_size_; }
    /// Returns pointer to the start of values for `node`-th node.
    template <typename T> T* access(std::vector<T>& v, int node) const {
        return v.data() + node*support_size_; }
    /// Returns const pointer to the start of values for `node`-th node for `op`-th operator.
    template <typename T> const T* access(const std::vector<T>& v, int op, int node) const {
        return v.data() + op*domain_size_*support_size_ + node*support_size_; }
    /// Returns const pointer to the start of values for `node`-th node.
    template <typename T> const T* access(const std::vector<T>& v, int node) const {
        return v.data() + node*support_size_; }
};

}  // namespace mm

#endif  // MEDUSA_BITS_OPERATORS_UNIFORMSHAPESTORAGE_FWD_HPP_
