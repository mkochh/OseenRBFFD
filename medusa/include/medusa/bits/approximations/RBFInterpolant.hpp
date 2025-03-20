#ifndef MEDUSA_BITS_APPROXIMATIONS_RBFINTERPOLANT_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_RBFINTERPOLANT_HPP_

/**
 * @file
 * Implementation of the RBFInterpolant class.
 */

#include "RBFInterpolant_fwd.hpp"
#include "RBFBasis.hpp"
#include <medusa/bits/utils/assert.hpp>

namespace mm {

template <typename RBFType, typename vec_t>
RBFInterpolant<RBFType, vec_t>::RBFInterpolant(
        const rbf_t& rbf, const Monomials<vec_t>& mon, const vector_t& point,
        const std::vector<vector_t>& support, scalar_t scale,
        const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>& coefficients) :
            rbf_(rbf), mon_(mon), point_(point), support_(support), scale_(scale),
            coefficients_(coefficients) {
    int n = support_.size();
    for (int i = 0; i < n; ++i) {
        support_[i] -= point;
        support_[i] /= scale;
    }
}

template <typename RBFType, typename vec_t>
typename vec_t::scalar_t RBFInterpolant<RBFType, vec_t>::operator()(const vector_t& point) const {
    int n1 = support_.size();
    int n2 = mon_.size();
    vector_t local_point = (point - point_) / scale_;
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> b(n1 + n2);
    for (int i = 0; i < n1; ++i) {
        b(i) = rbf_((local_point - support_[i]).squaredNorm());
    }
    for (int i = 0; i < n2; ++i) {
        b(n1+i) = mon_.eval(i, local_point);
    }
    return b.dot(coefficients_);
}

/// @cond
template <typename RBFType, typename vec_t>
template <typename operator_t>
typename vec_t::scalar_t
RBFInterpolant<RBFType, vec_t>::operator()(const vector_t& point, const operator_t& op) const {
    int n1 = support_.size();
    int n2 = mon_.size();
    vector_t local_point = (point - point_) / scale_;
    RBFBasis<rbf_t, vec_t> rbf_basis(n1, rbf_);
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> b(n1 + n2);
    for (int i = 0; i < n1; ++i) {
        b(i) = rbf_basis.evalOp(i, local_point, op, support_);
    }
    for (int i = 0; i < n2; ++i) {
        b(n1+i) = mon_.evalOp(i, local_point, op);
    }
    return b.dot(coefficients_);
}
/// @endcond

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_RBFINTERPOLANT_HPP_
