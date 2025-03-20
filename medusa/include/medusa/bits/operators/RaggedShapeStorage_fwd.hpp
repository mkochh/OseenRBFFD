#ifndef MEDUSA_BITS_OPERATORS_RAGGEDSHAPESTORAGE_FWD_HPP_
#define MEDUSA_BITS_OPERATORS_RAGGEDSHAPESTORAGE_FWD_HPP_

/**
 * @file
 * Declaration of shape storage for shapes of different sizes.
 *
 * @example test/operators/RaggedShapeStorage_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/operators/ShapeStorage_fwd.hpp>
#include <medusa/bits/operators/shape_flags.hpp>
#include <Eigen/Core>

namespace mm {

/**
 * Efficiently stores shape functions of different lengths.
 * This class is used to store shape functions (stencil weights)
 * generated for discretizations where supports of nodes have different sizes
 * e.g.\ using FindBalancedSupport. This class
 * is more efficient than storing the shapes in a nested type, such as
 * `std::vector<std::vector<T>>`, see
 * [technical report](http://http://e6.ijs.si/medusa/wiki/index.php/File:Tech_report.pdf).
 *
 * If supports of all nodes have the same size, use UniformShapeStorage instead.
 *
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
 * @snippet operators/RaggedShapeStorage_test.cpp Ragged shape storage usage example
 * @ingroup operators
 *
 * @sa sh, UniformShapeStorage, ShapeStorage, computeShapes
 */
template <typename vec_t, typename OpFamilies =
        std::tuple<Lap<vec_t::dim>, Der1s<vec_t::dim>, Der2s<vec_t::dim>>>
class RaggedShapeStorage :
        public ShapeStorage<RaggedShapeStorage<vec_t, OpFamilies>, vec_t, OpFamilies> {
  public:
    typedef vec_t vector_t;  ///< Vector type used.
    typedef typename vec_t::scalar_t scalar_t;  ///< Scalar type used.
    /// Store dimension of the domain.
    enum { /** Dimensionality of the domain. */ dim = vec_t::dim };

  private:
    /// Support sizes.
    Range<int> support_sizes_;
    /// Indexes of starts of supports. Cumulative sums of `support_sizes_`.
    Range<int> support_starts_;
    /// Sum of all support sizes.
    int total_size_;

    /// Parent class
    typedef ShapeStorage<RaggedShapeStorage<vec_t, OpFamilies>, vec_t, OpFamilies> base_t;
    friend base_t;  ///< Be friends with derived class.
    using base_t::domain_size_;
    using base_t::support_;
    using base_t::shapes_;

  public:
    using base_t::num_operators;  ///< Number of operators stored in this storage.
    using base_t::size;  ///< Number of nodes that shapes can be stored for.

    /// Constructs an empty shape storage with @ref size 0.
    RaggedShapeStorage() : total_size_(0) {}

    /// Returns support size of `node`-th node.
    int supportSize(int node) const { return support_sizes_[node]; }

  private:
    /**
     * Resizes the storage to accommodate shapes of given sizes. If support sizes are `{9, 12, 7}`
     * the class will allocate space for shapes for 3 nodes with sizes 9, 12 and 7.
     * The containers are zero initialized.
     */
    void resize_(const std::vector<int>& support_sizes);

    /// Returns pointer to the start of values for `node`-th node for `op`-th operator.
    template <typename T> T* access(std::vector<T>& v, int i, int j) const {
        return v.data() + i*total_size_ + support_starts_[j]; }
    /// Returns pointer to the start of values for `node`-th node for.
    template <typename T> T* access(std::vector<T>& v, int i) const {
        return v.data() + support_starts_[i]; }
    /// Returns const pointer to the start of values for `node`-th node for `op`-th operator.
    template <typename T> const T* access(const std::vector<T>& v, int i, int j) const {
        return v.data() + i*total_size_ + support_starts_[j]; }
    /// Returns const pointer to the start of values for `node`-th node.
    template <typename T> const T* access(const std::vector<T>& v, int i) const {
        return v.data() + support_starts_[i]; }
};

}  // namespace mm

#endif  // MEDUSA_BITS_OPERATORS_RAGGEDSHAPESTORAGE_FWD_HPP_
