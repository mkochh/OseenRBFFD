#ifndef MEDUSA_BITS_APPROXIMATIONS_RBFFD_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_RBFFD_HPP_

/**
 * @file
 * Implementation of Radial Basis Function Finite Difference approximation.
 */

#include "RBFFD_fwd.hpp"
#include <vector>
#include <medusa/bits/utils/numutils.hpp>
#include <Eigen/LU>
#include "RBFInterpolant.hpp"

namespace mm {

/// @cond
template <class RBFType, class vec_t, class scale_t, class solver_t>
Eigen::Matrix<typename vec_t::scalar_t, Eigen::Dynamic, Eigen::Dynamic>
RBFFD<RBFType, vec_t, scale_t, solver_t>::getMatrix(
        const std::vector<vector_t>& local_support) const {
    int n1 = local_support.size();
    assert_msg(n1 > 0, "Cannot construct a RBFFD approximation with no stencil nodes. "
                       "Did you forget to cal findSupport?.");
    int n2 = mon_.size();
    int N = n1 + n2;
    Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> M(N, N);
    // RBF part
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < i; ++j) {
            M(i, j) = M(j, i) = rbf_((local_support[i]-local_support[j]).squaredNorm());
        }
        M(i, i) = rbf_(0);
    }

    // Monomial part
    for (int i = 0; i < n2; ++i) {
        for (int j = 0; j < n1; ++j) {
            M(n1+i, j) = M(j, n1+i) = mon_.eval(i, local_support[j], local_support);
        }
    }
    M.bottomRightCorner(n2, n2).setZero();
    return M;
}

template <class RBFType, class vec_t, class scale_t, class solver_t>
void RBFFD<RBFType, vec_t, scale_t, solver_t>::compute(
        const vec_t& point, const std::vector<vec_t>& support) {
    scale_ = scale_t::scale(point, support);
    int n1 = support.size();
    // Local scaled support
    support_.resize(n1);
    for (int i = 0; i < n1; ++i) {
        support_[i] = (support[i] - point) / scale_;
    }
    Eigen::Matrix<typename vec_t::scalar_t, Eigen::Dynamic, Eigen::Dynamic> M = getMatrix(support_);
    point_ = point;
    solver_.compute(M);
}

template <class RBFType, class vec_t, class scale_t, class solver_t>
Eigen::Matrix<typename vec_t::scalar_t, Eigen::Dynamic, 1>
RBFFD<RBFType, vec_t, scale_t, solver_t>::getShape() const {
    int n = support_.size();
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> rhs(n + mon_.size());
    for (int i = 0; i < n; ++i) {
        rhs(i) = rbf_(support_[i].squaredNorm());
    }
    for (int i = 0; i < mon_.size(); ++i) {
        rhs(n+i) = mon_.evalAt0(i);
    }
    return getShapeFromRhs(rhs);
}

template <class RBFType, class vec_t, class scale_t, class solver_t>
template <class operator_t>
Eigen::Matrix<typename vec_t::scalar_t, Eigen::Dynamic, 1>
RBFFD<RBFType, vec_t, scale_t, solver_t>::getShape(const operator_t& op) const {
    int n = support_.size();
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> rhs(n + mon_.size());
    RBFBasis<RBFType, vec_t> rbf_basis(n, rbf_);
    for (int i = 0; i < n; ++i) {
        rhs(i) = rbf_basis.evalOpAt0(i, op, support_, scale_);
    }
    for (int i = 0; i < mon_.size(); ++i) {
        rhs(n+i) = mon_.evalOpAt0(i, op, support_, scale_);
    }
    return getShapeFromRhs(rhs);
}

template <class RBFType, class vec_t, class scale_t, class solver_t>
template <typename Derived>
Eigen::Matrix<typename vec_t::scalar_t, Eigen::Dynamic, 1>
RBFFD<RBFType, vec_t, scale_t, solver_t>::getShapeFromRhs(
        const Eigen::MatrixBase<Derived>& Lb) const {
    return solver_.solve(Lb).head(support_.size());
}

template <class RBFType, class vec_t, class scale_t, class solver_t>
template <typename Derived>
RBFInterpolant<RBFType, vec_t>
RBFFD<RBFType, vec_t, scale_t, solver_t>::getApproximant(const vector_t& point,
        const std::vector<vector_t>& support, const Eigen::MatrixBase<Derived>& values) const {
    int n = support.size();
    assert_msg(values.size() == n, "Number of given values %d must be equal to support size %d.",
               values.size(), n);
    double scale = scale_t::scale(point, support);
    // Local scaled support
    std::vector<vec_t> local_support(n);
    for (int i = 0; i < n; ++i) {
        local_support[i] = (support[i] - point) / scale;
    }
    Eigen::Matrix<typename vec_t::scalar_t, Eigen::Dynamic, Eigen::Dynamic> M =
        getMatrix(local_support);
    solver_t solver;
    solver.compute(M);
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> rhs(n + mon_.size());
    rhs.head(n) = values;
    rhs.tail(mon_.size()).setZero();
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> coeff = solver.solve(rhs);
    return RBFInterpolant<rbf_t, vec_t>(rbf_, mon_, point, support, scale, coeff);
}
/// @endcond

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_RBFFD_HPP_
