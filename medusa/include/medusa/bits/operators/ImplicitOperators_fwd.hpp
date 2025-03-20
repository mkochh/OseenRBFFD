#ifndef MEDUSA_BITS_OPERATORS_IMPLICITOPERATORS_FWD_HPP_
#define MEDUSA_BITS_OPERATORS_IMPLICITOPERATORS_FWD_HPP_

/**
 * @file
 * Declarations of implicit operators.
 *
 * @example test/operators/ImplicitOperators_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/approximations/Operators_fwd.hpp>
#include <medusa/bits/utils/assert.hpp>

namespace mm {

/**
 * This class represents implicit operators that fill given matrix @ref M and right hand side
 * @ref rhs with appropriate coefficients approximating differential operators with shape
 * functions from given shape storage @ref ss.
 *
 * @warning The operators grad(), value() etc.\ are lazy and will not write to the matrix
 * until summed with another operator or assigned to. If eager evaluation is required,
 * use the @ref OpBase::eval "eval()" method on each operation.
 *
 * @warning Each evaluated operation is only <b>added</b> to the matrix and rhs. Setting the same
 * equation twice thus has the effect of multiplying it by 2. Therefore, the matrix and the rhs
 * should be initialized (usually to zero) before setting any equations.
 *
 * @tparam shape_storage_type Any shape storage class satisfying the @ref ss-concept.
 * @tparam matrix_type Usually and Eigen sparse matrix.
 * @tparam rhs_type Usually an Eigen `VectorXd`.
 * @sa ImplicitVectorOperators
 *
 * Usage example:
 * @snippet operators/ImplicitOperators_test.cpp Implicit operators usage example
 * @ingroup operators
 */
template <class shape_storage_type, class matrix_type, class rhs_type>
class ImplicitOperators {
  public:
    typedef shape_storage_type shape_storage_t;  ///< Type of shape storage.
    typedef matrix_type matrix_t;  ///< Matrix type.
    typedef rhs_type rhs_t;  ///< Right hand side type.
    typedef typename shape_storage_t::vector_t vector_t;  ///< Vector type.
    /// Scalar type of matrix elements. Could be different from shape storage scalar type, e.g.
    /// shape storage type can be `double`, but the matrix could be `complex<double>`.
    typedef typename matrix_t::Scalar scalar_t;
    /// Store dimension of the domain.
    enum { /** Dimensionality of the function domain. */ dim = shape_storage_t::dim };

  private:
    /// Shape storage, but name is shortened for readability.
    const shape_storage_t* ss;
    matrix_t* M;  ///< Pointer to problem matrix.
    rhs_t* rhs;  ///< Pointer to right hand side.
    int row_offset;  ///< Row offset to be used when accessing matrix or rhs coefficients.
    int col_offset;  ///< Column offset to be used when accessing matrix coefficients.

    /**
     * Add `v` to `M(i, j)`.
     * @throw Assertion fails if `v` is none or `i` or `j` are out of range.
     */
    void addToM(int i, int j, scalar_t v);
    /**
     * Adds `v` to `rhs[i]`.
     * @throw Assertion fails if `v` is none or `i` is out of range.
     */
    void addToRhs(int i, scalar_t v);

    template <typename derived_t> class OpBase;  // Forward declaration needed.

    /**
     * Class representing an operation on a specific row of the matrix.
     * This row operations can be summed (with no effect) and assigned to, which sets
     * the right hand side. They enable a syntax like
     * `op.lap(i) + 2*op.grad(i, {1, 2}) = 2.3`.
     */
    class RowOp {
        friend class ImplicitOperators;
      protected:
        ImplicitOperators& op;  ///< Reference to underlying operators.
        int row_;  ///< Index of the row to which this operation refers. Row offset is not included.
        /// Construct a row operation for a given row.
        RowOp(ImplicitOperators& op, int row) : op(op), row_(row) {}
        /// Return absolute index of this row operation.
        int row() const { return op.row_offset + row_; }
      public:
        /**
         * Add two row operations together. This does nothing except for checking
         * that both operations refer to the same actual matrix row, catching some programming
         * errors.
         * @throws Assertion fails if `*this` and `other` do not refer to the same matrix row.
         */
        RowOp operator+(const RowOp& other) const {
            assert_msg(row() == other.row(),
                       "Cannot add together terms for different matrix rows %d and %d.",
                       row(), other.row());
            return *this;
        }
        /// Add normal operation to a row operation by evaluating the former and adding normally.
        template <typename derived_t>
        RowOp operator+(const OpBase<derived_t>& right) { return *this + right.eval(); }

        /// Assigns given value to the right hand side.
        void operator=(scalar_t value) {
            op.addToRhs(row_, value);
        }
    };

    /**
     * Base class for all elementary implicit operations. This class enables multiplication with
     * scalar, addition and evaluation of implicit operations. The method eval() is the
     * implementation of the actual operation that writes coefficients to the matrix.
     * CRTP is used for compile time polymorphism of `eval()`, so that concrete operations
     * only need to define this method.
     */
    template <typename derived_t>
    class OpBase {
      protected:
        ImplicitOperators& op;  ///< Reference to underlying operators.
        int node;  ///< Index of the point for which to apply the operation.
        int row;  ///< Matrix row to which the operation should be applied (without offset).
        /// Construct operation for a given node.
        OpBase(ImplicitOperators& op, int node, int row)
                : op(op), node(node), row(row) {}
      public:
        /**
         * Write appropriate coefficients for this operation to the matrix.
         * @param alpha Scalar to multiply the coefficients with.
         * @return A RowOp for the appropriate matrix row.
         */
        RowOp eval(scalar_t alpha) const {
            return static_cast<const derived_t &>(*this).eval(alpha);
        }
        /// Eval with `alpha = 1.0`.
        RowOp eval() const { return eval(1.0); }
        /// Evaluates `*this` and sets rhs to `x`.
        void operator=(scalar_t x) { eval() = x; }
        /// Multiply this operation by a scalar.
        RowOp operator*(scalar_t alpha) { return eval(alpha); }
        /// Combine this operation with another row operation.
        RowOp operator+(RowOp right) { return eval() + right; }
        /// Multiply this operation by `-1`.
        RowOp operator-() { return eval(-1.0); }
        /// Combine tho operation with another operation. Both operations are evaluated.
        template <typename other_derived_t>
        RowOp operator+(const OpBase<other_derived_t>& right) { return eval() + right.eval(); }
        /// Multiply with scalar from the left.
        friend RowOp operator*(scalar_t alpha, const OpBase& o) { return o.eval(alpha); }
    };

/// Macro that inherits appropriate members from parent OpBase class.
#define USE_BASE(Name) \
  private: \
    using OpBase<Name>::op; \
    using OpBase<Name>::node; \
    using OpBase<Name>::row; \
    using OpBase<Name>::OpBase; \
    friend class ImplicitOperators; \
  public: \
    using OpBase<Name>::eval; \
    using OpBase<Name>::operator=;

    /// Class representing the "evaluate" operation, i.e.\ the zero-th derivative.
    class ValueOp : public OpBase<ValueOp> {
        USE_BASE(ValueOp)
        /// Evaluate this operator.
        RowOp eval(scalar_t alpha) const {
            op.addToM(row, node, alpha);
            return RowOp(op, row);
        }
    };
    /// Class representing Basic operators i.e. Der1, Der2, Lap or user defined custom operators.
    template <typename op_family_t>
    class BasicOp : public OpBase<BasicOp<op_family_t>> {
        USE_BASE(BasicOp<op_family_t>)
      private:
        typename op_family_t::operator_t o;  ///< Basic operator with precomputed shape functions.
      public:
        /// Construct given basic operator 'o'.
        BasicOp(ImplicitOperators& op, int node, int row, typename op_family_t::operator_t o) :
                OpBase<BasicOp<op_family_t>>(op, node, row), o(o) {}
        /// Evaluate this operator.
        RowOp eval(scalar_t alpha) const {
            // There is some room for minor optimization here: the index does not need to be
            // computed for 1-element families. A special Op could be added for that,
            // but is probably not worth the code complexity.
            int idx = op_family_t::index(o);
            for (int i = 0; i < op.ss->supportSize(node); ++i) {
                op.addToM(row, op.ss->support(node, i), alpha *
                    op.ss->template get<op_family_t>(idx, node, i));
            }
            return RowOp(op, row);
        }
    };
    /**
     * Class representing the directional derivative (gradient) operation.
     * If @ref vec is denoted by @f$\vec v@f$, then this class represents
     * @f$(\vec v \cdot \nabla)@f$ operator evaluated at node with index @ref node.
     */
    class GradOp : public OpBase<GradOp> {
        USE_BASE(GradOp)
        vector_t vec;  ///< Vector representing the direction of differentiation.

        /// Construct directional derivative operator with given direction.
        GradOp(ImplicitOperators& op, int node, int row, vector_t vec) :
                OpBase<GradOp>(op, node, row), vec(vec) {}
      public:
        /// Evaluate this operator.
        RowOp eval(scalar_t alpha) const {
            for (int i = 0; i < op.ss->supportSize(node); ++i) {
                scalar_t value = 0;
                for (int var = 0; var < op.dim; ++var) {
                    value += vec[var] * op.ss->d1(var, node, i);
                }
                op.addToM(row, op.ss->support(node, i), alpha * value);
            }
            return RowOp(op, row);
        }
    };

  public:
    /// Default constructor sets offset to `0` and pointers to `nullptr`.
    ImplicitOperators() : ss(nullptr), M(nullptr), rhs(nullptr), row_offset(0), col_offset(0) {}
    /// Set only shape storage, the rest as default constructor.
    explicit ImplicitOperators(const shape_storage_t& ss) : ss(&ss), M(nullptr), rhs(nullptr),
                                                            row_offset(0), col_offset(0) {}
    /**
     * Set shape storage and problem matrix and rhs.
     * @param ss Class storing all computed shapes.
     * @param M Problem matrix. Must have at least `ss->size()` rows and cols.
     * @param rhs Problem right hand side. Must have at least `ss->size()` rows.
     * @param row_offset Instead of counting rows from 0, count them from row `row_offset`.
     * @param col_offset Instead of counting columns from 0, count them from row `col_offset`.
     *
     * @warning Matrices `M` and `rhs` have values only added to them and should be zero
     * initialized by the user. Since this is a common mistake, a warning is printed if this is not
     * the case when in debug mode.
     *
     * This class is usually constructed directly from shape storage using the
     * `implicitOperators()` member function.
     */
    ImplicitOperators(const shape_storage_t& ss, matrix_t& M, rhs_t& rhs,
                      int row_offset = 0, int col_offset = 0);

    /// Sets current matrix and right hand side.
    void setProblem(matrix_t& M, rhs_t& rhs, int row_offset = 0, int col_offset = 0) {
        ImplicitOperators::M = &M;
        ImplicitOperators::rhs = &rhs;
        setRowOffset(row_offset);
        setColOffset(col_offset);
    }
    /// Sets row offset for given matrix, treating it as is the first row had index `row_offset`.
    void setRowOffset(int row_offset) {
        assert_msg(0 <= row_offset, "Row offset cannot be negative, got %d.", row_offset);
        ImplicitOperators::row_offset = row_offset;
    }
    /// Sets col offset for given matrix, treating it as is the first column had index `col_offset`.
    void setColOffset(int col_offset) {
        assert_msg(0 <= col_offset, "Col offset cannot be negative, got %d.", col_offset);
        ImplicitOperators::col_offset = col_offset;
    }

    /// Returns `true` if operators have a non-null pointer to storage.
    bool hasShapes() const { return ss != nullptr; }
    /// Returns `true` if operators have a non-null pointer to problem matrix.
    bool hasMatrix() const { return M != nullptr; }
    /// Returns `true` if operators have a non-null pointer to problem right hand side.
    bool hasRhs() const { return rhs != nullptr; }

    /**
     * Sets implicit equation that value of a field is equal to some other value.
     * The code `alpha*op.value(node) = v` fills the `node`-th row of matrix `M`
     * so that `node`-th row of the equation
     * @f$ M u = v @f$
     * is a good approximation of the equation
     * @f[ \alpha u(p) = v(p), @f]
     * where @f$p@f$ is the `node`-th point.
     *
     * @param node Index of a node from 0 to `N` for which to write the equation for.
     * This means only the `node`-th row of the matrix is changed.
     * User must make sure that the matrix is large enough.
     * @param row Write equation in this specific row. Row with index `node` is chosen by default.
     * @return A class representing this operation for lazy evaluation purposes.
     */
    ValueOp value(int node, int row) { return ValueOp(*this, node, row); }
    /// Same as @ref value with row equal to current node.
    ValueOp value(int node) { return value(node, node); }

    /**
     * Sets implicit equation that Laplacian of a field @f$ u@f$ at `node`-th point
     * is equal to some value. The code `alpha*op.lap(node) = v` fills the `node`-th row of
     * matrix `M` so that `node`-th row of the equation
     * @f$ M u = v @f$
     * is a good approximation of the equation
     * @f[ \alpha \nabla^2 u(p) = v(p), @f]
     * where @f$p@f$ is the `node`-th point.
     *
     * @param node Index of a node from 0 to `N` for which to write the equation for.
     * This means only the `node`-th row of the matrix is changed.
     * User must make sure that the matrix is large enough.
     * @param row Write equation in this specific row. Row with index `node` is chosen by default.
     * @return A class representing this operation for lazy evaluation purposes.
     *
     * Example:
     * @snippet ImplicitOperators_test.cpp Implicit laplace 1d example
     */
    BasicOp<Lap<dim>> lap(int node, int row) { return BasicOp<Lap<dim>>(*this, node, row, {}); }
    /// Same as @ref lap with row equal to current node.
    BasicOp<Lap<dim>> lap(int node) { return lap(node, node); }

    /**
     * Sets implicit equation that gradient of a field @f$ u@f$ at `node`-th point
     * along `v` is equal to some value.
     * The code `alpha*op.grad(node, v) = r` fills the `node`-th row of
     * matrix `M` so that `node`-th row of the equation
     * @f$ M u = r @f$
     * is a good approximation of the equation
     * @f[ \alpha \vec{v} \cdot (\nabla u)(p) = r(p), @f]
     * where @f$p@f$ is the `node`-th point and @f$\vec{v}@f$ is the vector `v`.
     * Alternatively, this can be viewed as setting the
     * directional derivative along `v` to be equal to `r`.
     *
     * @param node Index of a node from 0 to `N` for which to write the equation for.
     * This means only `node`-th row of the matrix is changed.
     * User must make sure that the matrix is large enough.
     * @param v Vector to multiply the gradient with.
     * @param row Write equation in this specific row. Row with index `node` is chosen by default.
     * @return A class representing this operation for lazy evaluation purposes.
     *
     * Example:
     * @snippet ImplicitOperators_test.cpp Implicit Grad example
     */
    GradOp grad(int node, const vector_t& v, int row) { return GradOp(*this, node, row, v); }
    /// Same as @ref grad with row equal to current node.
    GradOp grad(int node, const vector_t& v) { return grad(node, v, node); }

    /**
     * Sets neumann boundary conditions in node `node`. The code
     * `alpha*op.neumann(node, normal) = v`
     * fills the `node`-th row of matrix `M` so that `node`-th row of the equation
     * @f$ M u = v @f$
     * is a good approximation of the equation
     * @f[ \alpha \dpar{u}{\vec{n}}(p) = v(p), @f]
     * where @f$p@f$ is the `node`-th point and @f$\vec{n}@f$ is the `unit_normal`.
     *
     * This is the same as using grad(), but has additional semantic meaning of setting the
     * boundary conditions.
     *
     * @return A class representing this operation for lazy evaluation purposes.
     *
     * @sa grad
     */
    GradOp neumann(int node, const vector_t& unit_normal, int row) {
        assert_msg(std::abs(unit_normal.norm() - 1.0) < 1e-13,
                   "Normal %s is not of unit length, got %.16f.", unit_normal, unit_normal.norm());
        return grad(node, unit_normal, row);
    }
    /// Same as @ref neumann with row equal to current node.
    GradOp neumann(int node, const vector_t& n) { return neumann(node, n, node); }

    /**
     * Sets implicit equation that derivative of a field @f$ u@f$ at `node`-th point
     * with respect to `var` is equal to some value.
     *
     * The code `alpha*op.der1(node, 0) = v`
     * fills the `node`-th row of matrix `M` so that `node`-th row of the equation
     * @f$ M u = v @f$
     * is a good approximation of the equation
     * @f[ \alpha \dpar{u}{x}(p) = v(p), @f]
     * where @f$p@f$ is the `node`-th point and @f$x@f$ is the `0`-th variable.
     *
     * @param node Index of a node from 0 to `N` for which to write the equation for.
     * This means only `node`-th row of the matrix is changed.
     * User must make sure that the matrix is large enough.
     * @param var Variable with respect to which to derive.
     * @param row Write equation in this specific row. Row with index `node` is chosen by default.
     *
     * @return A class representing this operation for lazy evaluation purposes.
     */
    BasicOp<Der1s<dim>> der1(int node, int var, int row) {
        return BasicOp<Der1s<dim>>(*this, node, row, {var});
    }
    /// Same as @ref der1 with row equal to current node.
    BasicOp<Der1s<dim>> der1(int node, int var) { return der1(node, var, node); }

    /**
     * Sets implicit equation that second derivative of a field @f$ u@f$ at `node`-th point
     * with respect to `varmin` and `varmax` is equal to some value.
     *
     * The code `alpha*op.der2(node, 0, 1) = v`
     * fills the `node`-th row of matrix `M` so that `node`-th row of the equation
     * @f$ M u = v @f$
     * is a good approximation of the equation
     * @f[ \alpha \dpar{^2u}{x\partial y}(p) = v(p), @f]
     * where @f$p@f$ is the `node`-th point and @f$x@f$ and @f$y@f$ are the `0`-th
     * and the `1`-st variables.
     *
     * @param node Index of a node from 0 to `N` for which to write the equation for.
     * This means only `node`-th row of the matrix is changed.
     * User must make sure that the matrix is large enough.
     * @param varmin Smaller of the two variables with respect to which to derive.
     * @param varmax Grater of the two variables with respect to which to derive.
     * @param row Write equation in this specific row. Row with index `node` is chosen by default.
     *
     * @return A class representing this operation for lazy evaluation purposes.
     */
    BasicOp<Der2s<dim>> der2(int node, int varmin, int varmax, int row) {
        return BasicOp<Der2s<dim>>(*this, node, row, {varmin, varmax});
    }
    /// Same as @ref der2 with row equal to current node.
    BasicOp<Der2s<dim>> der2(int node, int varmin, int varmax) {
        return der2(node, varmin, varmax, node);
    }


    /**
     * Add the weights for operator `o` to the appropriate elements in the matrix.
     * @param node
     * @param o Operator family
     * @param row Write in this matrix row.
     * @return This function is lazy and returns an "operation" which will write the
     * weights in the matrix, when evaluated.
     */
    template <typename op_family_t>
    BasicOp<op_family_t> apply(int node, typename op_family_t::operator_t o, int row) {
        return BasicOp<op_family_t>(*this, node, row, o);
    }

    /// Same as @ref apply with row equal to current node.
    template <typename op_family_t>
    BasicOp<op_family_t> apply(int node, typename op_family_t::operator_t o) {
        return apply<op_family_t>(node, o, node);
    }

    /// Overload for default-constructible operator.
    template <typename op_family_t>
    BasicOp<op_family_t> apply(int node) { return apply<op_family_t>(node, {}, node); }


    /// Output basic info about given operators.
    template <typename S, typename M, typename R>
    friend std::ostream& operator<<(std::ostream& os, const ImplicitOperators<S, M, R>& op);
};

}  // namespace mm

#endif  // MEDUSA_BITS_OPERATORS_IMPLICITOPERATORS_FWD_HPP_
