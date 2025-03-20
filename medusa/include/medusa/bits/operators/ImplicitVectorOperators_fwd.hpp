#ifndef MEDUSA_BITS_OPERATORS_IMPLICITVECTOROPERATORS_FWD_HPP_
#define MEDUSA_BITS_OPERATORS_IMPLICITVECTOROPERATORS_FWD_HPP_

/**
 * @file
 * Declarations of implicit vector operators.
 *
 * @example test/operators/ImplicitVectorOperators_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/utils/assert.hpp>
#include "ImplicitOperators_fwd.hpp"

namespace mm {

/**
 * This class represents implicit vector operators that fill given matrix @ref M and right hand side
 * @ref rhs with appropriate coefficients approximating differential operators with shape
 * functions from given shape storage @ref ss. They are intended for implicit solutions
 * of vector PDEs. Component specific equations can also be set using @ref eq.
 *
 * @warning The operators grad(), value() etc.\ are lazy and will not write to the matrix
 * until summed with another operator or assigned to. If eager evaluation is required,
 * use the @ref VecOpBase::eval "eval()" method on each operation.
 *
 * @warning Each evaluated operation is only <b>added</b> to the matrix and rhs. Setting the same
 * equation twice thus has the effect of multiplying it by 2. Therefore, the matrix and the rhs
 * should be initialized (usually to zero) before setting any equations.
 *
 * @tparam shape_storage_type Any shape storage class satisfying the @ref ss-concept.
 * @tparam matrix_type Usually and Eigen sparse matrix.
 * @tparam rhs_type Usually an Eigen `VectorXd`.
 *
 * The matrix and rhs should be of size at least `N*dim`, where `N` is the number of nodes in the
 * domain and `dim` its dimensionality.
 * @sa ImplicitOperators
 *
 * Usage example:
 * @snippet operators/ImplicitVectorOperators_test.cpp Implicit vector operators usage example
 * @ingroup operators
 */
template <class shape_storage_type, class matrix_type, class rhs_type>
class ImplicitVectorOperators {
  public:
    typedef shape_storage_type shape_storage_t;  ///< Type of shape storage.
    typedef matrix_type matrix_t;  ///< Matrix type.
    typedef rhs_type rhs_t;  ///< Right hand side type.
    /// Vector type for vectors in the domain, usually float or double.
    typedef typename shape_storage_t::vector_t domain_vector_t;
    /// Scalar type of matrix elements. Could be different from shape storage scalar type, e.g.
    /// shape storage type can be `double`, but the matrix could be `complex<double>`.
    typedef typename matrix_t::Scalar scalar_t;
    /// Store dimension of the domain.
    enum { /** Dimensionality of the domain. */ dim = shape_storage_t::dim };
    ///< Vector type for elements of vector fields, can be double of complex double.
    typedef Eigen::Matrix<scalar_t, dim, 1> vector_t;  ///< Vector type.

  private:
    /// Shape storage, but name is shortened for readability.
    const shape_storage_t* ss;
    matrix_t* M;  ///< Pointer to problem matrix.
    rhs_t* rhs;  ///< Pointer to right hand side.
    int row_offset;  ///< Row offset to be used when accessing matrix or rhs coefficients.
    int col_offset;  ///< Column offset to be used when accessing matrix coefficients.

    /// Adds `v` to `M(i, j)`, i.e.\ `M(i, j) += v`.
    void addToM(int i, int j, scalar_t v) {
        assert_msg(v == v, "You tried to set M(%d, %d) to NAN. Check your computes shapes.",
                   row_offset+i, col_offset+j);
        M->coeffRef(row_offset+i, col_offset+j) += v;
    }
    /// Adds components of `v` to  `i`, `i+N`, ..., `i+(dim-1)*N`-th rows of `rhs`.
    void addToRhs(int i, const vector_t& v) {
        assert_msg(!v.array().isNaN().any(),
                   "You tried to set rhs(%d) to NaN. Check you computed shapes.", row_offset+i);
        for (int d = 0; d < dim; ++d) {
            rhs->operator[](row_offset + d * ss->size() + i) += v[d];
        }
    }

    template <typename derived_t> class VecOpBase;  // Forward declaration needed.

    /**
     * Class representing an operation on a specific set of rows of the matrix. If `N` is the
     * number of nodes in the domain, then this class represents an operation on rows
     * `row`, `row+N`, ..., `row+(dim-1)*N`.
     * This row operations can be summed (with no effect) and assigned to, which sets
     * the right hand side. They enable a syntax like
     * `op.lap(i) + 2*op.grad(i, {1, 2}) = 2.3`.
     */
    class RowVecOp {
        friend class ImplicitVectorOperators;
      protected:
        ImplicitVectorOperators& op;  ///< Reference to underlying operators.
        int row_;  ///< Index of the row to which this operation refers. Row offset is not included.
        /// Construct a row operation for a given row.
        RowVecOp(ImplicitVectorOperators& op, int row) : op(op), row_(row) {}
        /// Return absolute index of this row operation.
        int row() const { return op.row_offset + row_; }
      public:
        /**
         * Add two row operations together. This does nothing except for checking
         * that both operations refer to the same actual matrix row, catching some programming
         * errors.
         * @throws Assertion fails if `*this` and `other` do not refer to the same matrix row.
         */
        RowVecOp operator+(const RowVecOp& other) const {
            assert_msg(row() == other.row(),
                       "Cannot add together terms for different matrix rows %d and %d.",
                       row(), other.row());
            return *this;
        }
        /// Add normal operation to a row operation by evaluating the former and adding normally.
        template <typename derived_t>
        RowVecOp operator+(const VecOpBase<derived_t>& right) { return *this + right.eval(); }
        /// Assigns given value to the right hand side.
        void operator=(const vector_t& value) {
            op.addToRhs(row_, value);
        }
    };

    /**
     * Base class for all elementary implicit vector operations. This class enables multiplication
     * with scalar, addition and evaluation of implicit operations. The method eval() is the
     * implementation of the actual operation that writes coefficients to the matrix.
     * CRTP is used for compile time polymorphism of eval(), so that concrete operations
     * only need to define this method.
     */
    template <typename derived_t>
    class VecOpBase {
      protected:
        ImplicitVectorOperators& op;  ///< Reference to underlying operators.
        int node;  ///< Index of the point for which to apply the operation.
        int row;  ///< Matrix row to which the operation should be applied (without offset).
        /// Construct operation for a given node.
        VecOpBase(ImplicitVectorOperators& op, int node, int row)
                : op(op), node(node), row(row) {}
      public:
        /**
         * Write appropriate coefficients for this operation to the matrix.
         * @param alpha Scalar to multiply the coefficients with.
         * @return A RowVecOp for the appropriate matrix row.
         */
        RowVecOp eval(scalar_t alpha) const {
            return static_cast<const derived_t &>(*this).eval(alpha);
        }
        /// Eval with `alpha = 1.0`.
        RowVecOp eval() const { return eval(1.0); }
        /// Evaluates `*this` and sets rhs to `x`.
        void operator=(const vector_t& x) { eval() = x; }
        /// Multiply this operation by a scalar.
        RowVecOp operator*(scalar_t alpha) { return eval(alpha); }
        /// Combine this operation with another row operation.
        RowVecOp operator+(RowVecOp right) { return eval() + right; }
        /// Multiply this operation by `-1`.
        RowVecOp operator-() { return eval(-1.0); }
        /// Combine tho operation with another operation. Both operations are evaluated.
        template <typename other_derived_t>
        RowVecOp operator+(const VecOpBase<other_derived_t>& right) { return eval()+right.eval(); }
        /// Multiply with scalar from the left.
        friend RowVecOp operator*(scalar_t alpha, const VecOpBase& op) { return op.eval(alpha); }
    };

/// Macro that inherits appropriate members from parent VecOpBase class.
#define USE_VECTOR_BASE(Name) \
  private: \
    using VecOpBase<Name>::op; \
    using VecOpBase<Name>::node; \
    using VecOpBase<Name>::row; \
    using VecOpBase<Name>::VecOpBase; \
    friend class ImplicitVectorOperators; \
  public: \
    using VecOpBase<Name>::eval; \
    using VecOpBase<Name>::operator=;

    /// Class representing the "evaluate" operation, i.e.\ the zero-th derivative.
    class ValueOp : public VecOpBase<ValueOp> {
        USE_VECTOR_BASE(ValueOp)
        /// Evaluate this operator.
        RowVecOp eval(scalar_t alpha) const {
            for (int d = 0; d < dim; ++d) {
                op.addToM(row + d*op.ss->size(), node + d*op.ss->size(), alpha);
            }
            return RowVecOp(op, row);
        }
    };

    /// Class representing basic operators.
    template <typename op_family_t>
    class BasicOp : public VecOpBase<BasicOp<op_family_t>> {
      USE_VECTOR_BASE(BasicOp<op_family_t>)
      private:
        typename op_family_t::operator_t o;  ///< Basic operator with precomputed shape functions.
      public:
        /// Construct given basic operator 'o'.
        BasicOp(ImplicitVectorOperators& op, int node, int row, typename op_family_t::operator_t o)
                : VecOpBase<BasicOp<op_family_t>>(op, node, row), o(o) {}
        /// Evaluate this operator.
        RowVecOp eval(scalar_t alpha) const {
            int idx = op_family_t::index(o);
            for (int d = 0; d < dim; ++d) {
                for (int i = 0; i < op.ss->supportSize(node); ++i) {
                    op.addToM(d*op.ss->size() + row, d*op.ss->size() + op.ss->support(node, i),
                      alpha * op.ss->template get<op_family_t>(idx, node, i));
                }
            }
            return RowVecOp(op, row);
        }
    };

    /**
     * Class representing the directional derivative (gradient) operation.
     * If @ref vec is denoted by @f$\vec v@f$, then this class represents
     * @f$(\vec v \cdot \nabla)@f$ operator evaluated at node with index @ref node.
     */
    class GradOp : public VecOpBase<GradOp> {
        USE_VECTOR_BASE(GradOp)
      private:
        domain_vector_t vec;  ///< Vector representing the direction of differentiation.
        /// Construct directional derivative operator with given direction.
        GradOp(ImplicitVectorOperators& op, int node, int row, const domain_vector_t& vec) :
                VecOpBase<GradOp>(op, node, row), vec(vec) {}
      public:
        /// Evaluate this operator.
        RowVecOp eval(scalar_t alpha) const {
            for (int i = 0; i < op.ss->supportSize(node); ++i) {
                scalar_t shape_value = 0;
                for (int var = 0; var < dim; ++var) {
                    shape_value += vec[var] * op.ss->d1(var, node, i);
                }
                for (int d = 0; d < dim; ++d) {
                    op.addToM(d*op.ss->size() + row, d*op.ss->size() + op.ss->support(node, i),
                              alpha*shape_value);
                }
            }
            return RowVecOp(op, row);
        }
    };

    /**
     * Class representing the gradient of the divergence operator, i.e.\ operator
     * @f$\nabla\nabla\cdot@f$.
     */
    class GradDivOp : public VecOpBase<GradDivOp> {
        USE_VECTOR_BASE(GradDivOp)
        /// Evaluate this operator.
        RowVecOp eval(scalar_t alpha) const {
            // Implicit equation is for all d: \sum_k \sum_j \chi_djk u_j(s_k) = v_d(s)
            for (int d = 0; d < dim; ++d) {
                for (int k = 0; k < op.ss->supportSize(node); ++k) {
                    for (int j = 0; j < dim; ++j) {  // loop over dimensions
                        int dmin = std::min(d, j);
                        int dmax = std::max(d, j);
                        op.addToM(d*op.ss->size() + row,
                                  j*op.ss->size() + op.ss->support(node, k),
                                  alpha * op.ss->d2(dmin, dmax, node, k));
                    }
                }
            }
            return RowVecOp(op, row);
        }
    };

    /// Class representing the traction operator. Useful for setting boundary conditions.
    class TractionOp : public VecOpBase<TractionOp> {
      USE_VECTOR_BASE(TractionOp)
      private:
        scalar_t lam,  ///<  1st Lame constant.
                 mu;  ///< 2nd Lame constant.
        domain_vector_t normal;  ///< Unit normal to the surface where traction os desired.

        /// Create traction operator with given Lame constants and outside unit normal.
        TractionOp(ImplicitVectorOperators& op, int node, int row, scalar_t lam, scalar_t mu,
                   const domain_vector_t& normal) :
                VecOpBase<TractionOp>(op, node, row), lam(lam), mu(mu), normal(normal) {}
      public:
        /// Evaluate this operator.
        RowVecOp eval(scalar_t alpha) const {
            for (int i = 0; i < dim; ++i) {  // which equation of sigma.n = t
                for (int j = 0; j < dim; ++j) {  // computation of sigma (lam tr(eps) I) part
                    op.eq(i).c(j).der1(node, j, row).eval(alpha*lam*normal(i));
                }
                for (int j = 0; j < dim; ++j) {  // computation of sigma (2 mu eps) part
                    op.eq(i).c(i).der1(node, j, row).eval(alpha*mu*normal(j));
                    op.eq(i).c(j).der1(node, i, row).eval(alpha*mu*normal(j));
                }
            }
            return RowVecOp(op, row);
        }
    };

    /// Represents one scalar component of a vector equation.
    class Equation {
        ImplicitVectorOperators& op;  ///< Reference to underlying operators.
        int eq_ixd;  ///< Index of the equation.
        friend class ImplicitVectorOperators;

        /// Create an equation with given index.
        Equation(ImplicitVectorOperators& op, int idx) : op(op), eq_ixd(idx) {
            assert_msg(0 <= idx && idx < op.dim, "Equation number %d out of range [0, %d).",
                       idx, dim);
        }

      public:
        /// Returns ordinary ImplicitOperators for component `comp` of equation `eq_idx`. @sa eq
        ImplicitOperators<shape_storage_t, matrix_t, rhs_t> c(int comp) {
            assert_msg(0 <= comp && comp < dim,
                       "Dimension %d out of range, expected in range [0, %d).", comp, dim);
            auto sop = ImplicitOperators<shape_storage_t, matrix_t, rhs_t>(*op.ss);
            sop.setProblem(*op.M, *op.rhs, op.row_offset + eq_ixd*op.ss->size(),
                    op.col_offset + comp*op.ss->size());
            return sop;
        }
    };

  public:
    /// Default constructor sets offset to `0` and pointers to `nullptr`.
    ImplicitVectorOperators() : ss(nullptr), M(nullptr), rhs(nullptr), row_offset(0),
                                col_offset(0) {}
    /// Set only shape storage, the rest as default constructor.
    explicit ImplicitVectorOperators(const shape_storage_t& ss) : ss(&ss), M(nullptr), rhs(nullptr),
                                                                  row_offset(0), col_offset(0) {}
    /**
     * Set shape storage and problem matrix and rhs.
     * @param ss Class storing all computed shapes.
     * @param M Problem matrix. Must have at least `ss->size()*::dim` rows and cols.
     * @param rhs Problem right hand side. Must have at least `ss->size()*::dim` rows.
     * @param row_offset Instead of counting rows from 0, count them from row `row_offset`.
     * @param col_offset Instead of counting columns from 0, count them from row `col_offset`.
     *
     * @warning Matrices `M` and `rhs` have values only added to them and should be zero
     * initialized by the user. Since this is a common mistake, a warning is printed if this is not
     * the case when in debug mode.
     *
     * This class is usually constructed directly from shape storage using the
     * `implicitVectorOperators()` member function.
     */
    ImplicitVectorOperators(const shape_storage_t& ss, matrix_t& M, rhs_t& rhs,
                            int row_offset = 0, int col_offset = 0);

    /// Sets current matrix and right hand side.
    void setProblem(matrix_t& M, rhs_t& rhs, int row_offset = 0, int col_offset = 0) {
        ImplicitVectorOperators::M = &M;
        ImplicitVectorOperators::rhs = &rhs;
        setRowOffset(row_offset);
        setColOffset(col_offset);
    }
    /// Sets row offset for given matrix, treating it as is the first row had index `row_offset`.
    void setRowOffset(int row_offset) {
        assert_msg(0 <= row_offset, "Row offset cannot be negative, got %d.", row_offset);
        ImplicitVectorOperators::row_offset = row_offset;
    }
    /// Sets col offset for given matrix, treating it as is the first column had index `col_offset`.
    void setColOffset(int col_offset) {
        assert_msg(0 <= col_offset, "Col offset cannot be negative, got %d.", col_offset);
        ImplicitVectorOperators::col_offset = col_offset;
    }

    /// Returns `true` if operators have a non-null pointer to storage and `false` otherwise.
    bool hasShapes() const { return ss != nullptr; }
    /// Returns `true` if operators have a non-null pointer to problem matrix.
    bool hasMatrix() const { return M != nullptr; }
    /// Returns `true` if operators have a non-null pointer to problem right hand side.
    bool hasRhs() const { return rhs != nullptr; }
    /**
     * Choose one specific component of the vector equation to write to.
     * Usage example:
     * @snippet ImplicitVectorOperators_test.cpp Eq usage example
     * sets equation @f$\alpha \dpar{u}{x}(i) + \beta \dpar{v}{y}(i) = t_x(i)@f$ as the first scalar
     * equation in the matrix.
     * @sa Equation::c
     */
    Equation eq(int num) {
        assert_msg(0 <= num && num < dim, "Equation must be in range [0, %d), got %d", num, dim);
        return Equation(*this, num);
    }

    /**
     * Sets implicit equation that value of a vector field is equal to some other value.
     *
     * The code `alpha*op.value(node) = v`
     * fills the `node`-th, `node+N`-th and `node+2*N`-th row of matrix `M`
     * these rows of equation @f$ M u = v @f$ are a good approximation of the equation
     * @f[ \alpha \vec{u}(p) = \vec{v}(p), @f]
     * where @f$p@f$ is the `node`-th point.
     *
     * @note Above example was made for `dim = 3`, other dimensions are analogous.
     * @param node Index of a node from 0 to `N` for which to write the equation for.
     * This means only `node`-th`, node+N`-th and `node+2*N`-th rows of the matrix are changed.
     * User must make sure that the matrix is large enough.
     * @param row Write equation in this specific row. Row with index `node` is chosen by default.
     *
     * @return A class representing this operation for lazy evaluation purposes.
     */
    ValueOp value(int node, int row) { return ValueOp(*this, node, row); }
    /// Same as @ref value with row equal to current node.
    ValueOp value(int node) { return value(node, node); }

    /**
     * Sets implicit equation that Laplacian of a vector field is equal to some other value.
     *
     * The code `alpha*op.lap(node) = v`
     * fills the `node`-th, `node+N`-th and `node+2*N`-th row of matrix `M`
     * these rows of equation @f$ M u = v @f$ are a good approximation of the equation
     * @f[ \alpha \nabla^2\vec{u}(p) = \vec{v}(p), @f]
     * where @f$p@f$ is the `node`-th point.
     *
     * @note Above example was made for `dim = 3`, other dimensions are analogous.
     * @param node Index of a node from 0 to `N` for which to write the equation for.
     * This means only `node`-th`, node+N`-th and `node+2*N`-th rows of the matrix are changed.
     * User must make sure that the matrix is large enough.
     * @param row Write equation in this specific row. Row with index `node` is chosen by default.
     *
     * @return A class representing this operation for lazy evaluation purposes.
     *
     * Example:
     * @snippet ImplicitVectorOperators_test.cpp Implicit vector laplace 2d example
     */
    BasicOp<Lap<dim>> lap(int node, int row) { return BasicOp<Lap<dim>>(*this, node, row); }
    /// Same as @ref lap with row equal to current node.
    BasicOp<Lap<dim>> lap(int node) { return lap(node, node); }

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

    /**
     * Sets implicit equation that gradient of a vector field along `v` is equal to some other
     * value. The code `alpha*op.grad(node, v) = r`
     * fills the `node`-th, `node+N`-th and `node+2*N`-th row of matrix `M`
     * these rows of equation @f$ M u = r @f$ are a good approximation of the equation
     * @f[ \alpha (\vec{v} \cdot \nabla) \vec{u}(p) = \vec{r}(p), @f]
     * where @f$p@f$ is the `node`-th point and @f$\vec{v}@f$ is the vector `v`.
     *
     * @note Above example was made for `dim = 3`, other dimensions are analogous.
     * @param node Index of a node from 0 to `N` for which to write the equation for.
     * This means only `node`-th`, node+N`-th and `node+2*N`-th rows of the matrix are changed.
     * User must make sure that the matrix is large enough.
     * @param row Write equation in this specific row. Row with index `node` is chosen by default.
     * @param v Vector to multiply the gradient with.
     *
     * Example:
     * @snippet ImplicitVectorOperators_test.cpp Implicit Gradvec example
     */
    GradOp grad(int node, const domain_vector_t& v, int row) { return GradOp(*this, node, row, v); }
    /// Same as @ref grad with row equal to current node.
    GradOp grad(int node, const domain_vector_t& v) { return grad(node, v, node); }

    /**
     * Sets implicit equation that gradient of divergence of a vector field is equal to some other
     * value. The code `alpha*op.graddiv(node) = v`
     * fills the `node`-th, `node+N`-th and `node+2*N`-th row of matrix `M`
     * these rows of equation @f$ M u = v @f$ are a good approximation of the equation
     * @f[ \alpha \nabla(\nabla \cdot \vec{u})(p) = \vec{v}(p), @f]
     * where @f$p@f$ is the `node`-th point.
     *
     * @note Above example was made for `dim = 3`, other dimensions are analogous.
     * @param node Index of a node from 0 to `N` for which to write the equation for.
     * This means only `node`-th`, node+N`-th and `node+2*N`-th rows of the matrix are changed.
     * User must make sure that the matrix is large enough.
     * @param row Write equation in this specific row. Row with index `node` is chosen by default.
     *
     * Example:
     * @snippet ImplicitVectorOperators_test.cpp Implicit graddiv example
     */
    GradDivOp graddiv(int node, int row) { return GradDivOp(*this, node, row); }
    /// Same as @ref graddiv with row equal to current node.
    GradDivOp graddiv(int node) { return graddiv(node, node); }

    /**
     * Sets neumann boundary conditions in node `node`. The code
     * `alpha*op.neumann(node, normal) = v`
     * fills the `node`-th, `node+N`-th and `node+2*N`-th row of matrix `M`
     * these rows of equation @f$ M u = v @f$ are a good approximation of the equation
     * @f[ \alpha \dpar{\vec{u}}{\vec{n}}(p) = \vec{v}(p) @f]
     * where @f$p@f$ is the `node`-th point and @f$\vec{n}@f$ is the `unit_normal`.
     *
     * This is the same as using grad(), but has additional semantic meaning of setting the
     * boundary conditions.
     *
     * @return A class representing this operation for lazy evaluation purposes.
     *
     * @sa grad
     */
    GradOp neumann(int node, const domain_vector_t& unit_normal, int row) {
        assert_msg(std::abs(unit_normal.norm() - 1.0) < 1e-13,
                   "Normal %s is not of unit length, got %.16f.", unit_normal, unit_normal.norm());
        return grad(node, unit_normal, row);
    }
    /// Same as @ref neumann with row equal to current node.
    GradOp neumann(int node, const domain_vector_t& n) { return neumann(node, n, node); }

    /**
     * Sets traction boundary conditions in node `node`. The code
     * `alpha*op.traction(node, lam, mu, normal) = t`
     * fills the `node`-th, `node+N`-th and `node+2*N`-th row of matrix `M`
     * these rows of equation @f$ M u = t @f$ are a good approximation of the equation
     * @f[ \alpha (\sigma(\vec{u}) \vec{n}) (p) = \vec{t}(p) @f]
     * where @f$p@f$ is the `node`-th point, @f$\sigma@f$ is the stress tensor
     * and @f$\vec{n}@f$ is the `unit_normal`.
     *
     * @param node Index of a node from 0 to `N` for which to write the equation for.
     * This means only `node`-th`, node+N`-th and `node+2*N`-th rows of the matrix are changed.
     * User must make sure that the matrix is large enough.
     * @param row Write equation in this specific row. Row with index `node` is chosen by default.
     * @param lam The first Lame coefficient.
     * @param mu The second Lame coefficient.
     * @param unit_normal The outside unit normal on the surface for which to set the traction.
     *
     * @return A class representing this operation for lazy evaluation purposes.
     */
    TractionOp traction(int node, scalar_t lam, scalar_t mu, const domain_vector_t& unit_normal,
                        int row) {
        assert_msg(std::abs(unit_normal.norm() - 1.0) < 1e-13,
                   "Normal %s is not of unit length, got %.16f.", unit_normal, unit_normal.norm());
        return TractionOp(*this, node, row, lam, mu, unit_normal);
    }
    /// Same as @ref traction with row equal to current node.
    TractionOp traction(int node, scalar_t lam, scalar_t mu, const domain_vector_t& n) {
        return traction(node, lam, mu, n, node);
    }

    /// Output basic info about given operators.
    template <typename S, typename M, typename R>
    friend std::ostream& operator<<(std::ostream& os, const ImplicitVectorOperators<S, M, R>& op);
};

}  // namespace mm

#endif  // MEDUSA_BITS_OPERATORS_IMPLICITVECTOROPERATORS_FWD_HPP_
