#ifndef MEDUSA_BITS_OPERATORS_SHAPESTORAGE_FWD_HPP_
#define MEDUSA_BITS_OPERATORS_SHAPESTORAGE_FWD_HPP_

/**
 * @file
 * Base class and interface for storing and accessing computed shapes.
 */

#include <medusa/Config.hpp>
#include <medusa/bits/types/Range_fwd.hpp>
#include <medusa/bits/utils/stdtypesutils.hpp>
#include <medusa/bits/approximations/Operators_fwd.hpp>
#include <iosfwd>
#include <tuple>
#include <Eigen/Core>

namespace mm {

template <class shape_storage_type>
class ExplicitOperators;

template <class shape_storage_type>
class ExplicitVectorOperators;

template <class shape_storage_type, class matrix_type, class rhs_type>
class ImplicitOperators;

template <class shape_storage_type, class matrix_type, class rhs_type>
class ImplicitVectorOperators;

/**
 * Shape storage base class. It supports storing computed shapes for Lap, D1 and D2 operators.
 * Derived classes should inherit the interface via CRTP and implement missing methods.
 *
 * @ingroup operators
 *
 * @sa UniformShapeStorage, RaggedShapeStorage
 */
template <typename Derived, typename vec_t, typename OpFamilies =
        std::tuple<Lap<vec_t::dim>, Der1s<vec_t::dim>, Der2s<vec_t::dim>>>
class ShapeStorage {
  public:
    typedef vec_t vector_t;  ///< Vector type used.
    typedef typename vec_t::scalar_t scalar_t;  ///< Scalar type used.
    /// Bitmask telling us which shapes to create.
    /// Store dimension of the domain.
    enum { /** Dimensionality of the domain. */ dim = vec_t::dim };
    typedef OpFamilies op_families_tuple;  ///< Tuple of operator families
    /// Number of operator families in this storage.
    constexpr static int num_operators = std::tuple_size<op_families_tuple>::value;

    friend Derived;  ///< Be friends with derived class.

  private:
    /// Total number of nodes.
    int domain_size_;
    /// Local copy of support domains.
    Range<int> support_;

    /// Tuple of shape containers for given operators.
    std::array<Range<scalar_t>, num_operators> shapes_;

    /// Returns pointer to the start of values for `node`-th node for `op`-th operator.
    template <typename T> SINL T* access(std::vector<T>& v, int op, int node) const {
        return static_cast<const Derived*>(this)->access(v, op, node); }
    /// Returns pointer to the start of values for `node`-th node.
    template <typename T> SINL T* access(std::vector<T>& v, int node) const {
        return static_cast<const Derived*>(this)->access(v, node); }
    /// Returns const pointer to the start of values for `node`-th node for `op`-th operator.
    template <typename T> SINL const T* access(const std::vector<T>& v, int op, int node) const {
        return static_cast<const Derived*>(this)->access(v, op, node); }
    /// Returns const pointer to the start of values for `node`-th node.
    template <typename T> SINL const T* access(const std::vector<T>& v, int node) const {
        return static_cast<const Derived*>(this)->access(v, node); }

    /// Constructs empty storage with size 0. Can be resized with @ref resize.
    ShapeStorage() : domain_size_(0) {}

  public:
    /// Resizes the storage to accommodate shapes of given sizes.
    void resize(const std::vector<int>& support_sizes) {
        domain_size_ = support_sizes.size();
        static_cast<Derived*>(this)->resize_(support_sizes);
    }

    /// Returns number of nodes.
    int size() const { return domain_size_; }

    /**
     * Returns index of `j`-th neighbour of `node`-th node, i.e.\ `support[node][j]`,
     * but possibly faster.
     */
    SINL int support(int node, int j) const;

    /// Returns a collection of support indices for given `node`.
    SINL Eigen::Map<const Eigen::Matrix<int, Eigen::Dynamic, 1>> support(int node) const;

    /// Sets support of `node`-th node to `support`.
    void setSupport(int node_idx, const std::vector<int>& support);

    /// Returns support size of `node`-th node.
    int supportSize(int node) const { return static_cast<const Derived*>(this)->supportSize(node); }

    /**
     * Returns a vector of support sizes for all nodes, useful for matrix space preallocation.
     * @warning If there are ghost nodes present, their support sizes might be zero. However,
     * matrix rows corresponding to ghost nodes must still be preallocated.
     * @sa supportSizesVec
     */
    Range<int> supportSizes() const;

    /**
     * Returns a `dim*N` vector of `dim*support_size`, useful for matrix space preallocation
     * for vector equations. Same warning as for @ref supportSizes applies.
     * @sa supportSizes
     */
    Range<int> supportSizesVec() const;

    // GENERAL GETTERS
    /// Read access to complete shape storage.
    const std::array<Range<scalar_t>, num_operators>& shapes() const { return shapes_; }

    /// Get the weight of `j`-th stencil node of `node` for operator `op` in `op_family` (by index).
    template <int op_family> SINL scalar_t get(int op, int node, int j) const {
        return access(shapes_[op_family], op, node)[j];
    }

    /// Get the weight of `j`-th stencil node of `node` for operator `op` in `op_family` (by type).
    template <typename op_family> SINL scalar_t get(int op, int node, int j) const {
        return get<tuple_index<op_family, op_families_tuple>::value>(op, node, j);
    }

    /// Get the `j`-th weight of shape for `node` for the only operator in `op_family` (by index).
    template <int op_family> SINL scalar_t get(int node, int j) const {
        return access(shapes_[op_family], node)[j];
    }

    /// Get the `j`-th weight of shape for `node` for the only operator in `op_family` (by type).
    template <typename op_family> SINL scalar_t get(int node, int j) const {
        return get<tuple_index<op_family, op_families_tuple>::value>(node, j);
    }

    /// Get the shape for `node` for operator `op` in `op_family` (by index).
    template <int op_family> Eigen::Map<const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>>
    getShape(int op, int node) const {
        return {access(shapes_[op_family], op, node), supportSize(node)};
    }

    /// Get the shape for `node` for operator `op` in `op_family` (by type).
    template <typename op_family> Eigen::Map<const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>>
    getShape(int op, int node) const {
        return getShape<tuple_index<op_family, op_families_tuple>::value>(op, node);
    }

    /// Get the shape for `node` for the only operator in `op_family` (by index).
    template <int op_family>
    SINL Eigen::Map<const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>> getShape(int node) const {
        return {access(shapes_[op_family], node), supportSize(node)};
    }

    /// Get the shape for `node` for the only operator in `op_family` (by type).
    template <typename op_family>
    SINL Eigen::Map<const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>> getShape(int node) const {
        return getShape<tuple_index<op_family, op_families_tuple>::value>(node);
    }

    // GENERAL SETTERS
    /// Set the shape for `node` for operator `op` in `op_family` (by index).
    template <int op_family>
    SINL void setShape(int op, int node, const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>& shape) {
        std::memcpy(access(shapes_[op_family], op, node), shape.data(),
                    supportSize(node)*sizeof(scalar_t));
    }

    /// Set the shape for `node` for operator `op` in `op_family` (by type).
    template <typename op_family>
    SINL void setShape(int op, int node, const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>& shape) {
        return setShape<tuple_index<op_family, op_families_tuple>::value>(op, node, shape);
    }

    /// Set the shape for `node` for the only operator in `op_family` (by index).
    template <int op_family>
    SINL void setShape(int node, const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>& shape) {
        std::memcpy(access(shapes_[op_family], node), shape.data(),
                    supportSize(node)*sizeof(scalar_t));    }

    /// Set the shape for `node` for the only operator in `op_family` (by type).
    template <typename op_family>
    SINL void setShape(int node, const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>& shape) {
        return setShape<tuple_index<op_family, op_families_tuple>::value>(node, shape);
    }

    // SPECIFIC ACCESSORS

    /// Return `j-th` laplace shape coefficient for `node`-th node.
    SINL scalar_t laplace(int node, int j) const;

    /// Returns the laplace shape for `node`.
    Eigen::Map<const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>> laplace(int node) const;

    /// Sets the laplace shape for `node` to `shape`.
    SINL void setLaplace(int node, const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>& shape);

    /// Return `j`-th shape coefficient for derivative wrt. variable `var` in `node`.
    SINL scalar_t d1(int var, int node, int j) const;

    /// Return shape for derivative wrt. variable `var` in `node`.
    Eigen::Map<const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>> d1(int var, int node) const;

    /// Sets shape for derivative wrt. variable `var` for `node` to `shape`.
    SINL void setD1(int var, int node, const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>& shape);

    /**
     * Return `j`-th shape coefficient for mixed derivative wrt. variables `varmin` and `varmax` in
     * `node`.
     */
    SINL scalar_t d2(int varmin, int varmax, int node, int j) const;

    /// Return shape for mixed derivative wrt. variables `varmin` and `varmax` in `node`.
    Eigen::Map<const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>>
    d2(int varmin, int varmax, int node) const;

    /// Sets shape for mixed derivative wrt. variables `varmin` and `varmax` for `node` to `shape`.
    SINL void setD2(int varmin, int varmax, int node,
                    const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>& shape);

    /// Returns the approximate memory used (in bytes).
    size_t memoryUsed() const;

    /// Construct explicit operators over this storage.
    ExplicitOperators<Derived> explicitOperators() const;
    /// Construct explicit vector operators over this storage.
    ExplicitVectorOperators<Derived> explicitVectorOperators() const;
    /// Construct implicit operators over this storage.
    template <typename M, typename R>
    ImplicitOperators<Derived, M, R> implicitOperators(M& matrix, R& rhs) const;
    /// Construct implicit vector operators over this storage.
    template <typename M, typename R>
    ImplicitVectorOperators<Derived, M, R> implicitVectorOperators(M& matrix, R& rhs) const;

    /// Output basic info about this shape storage.
    template <typename D, typename V, typename O>
    friend std::ostream& operator<<(std::ostream& os, const ShapeStorage<D, V, O>& shapes);
};


}  // namespace mm

#endif  // MEDUSA_BITS_OPERATORS_SHAPESTORAGE_FWD_HPP_
