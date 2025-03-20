#ifndef MEDUSA_BITS_APPROXIMATIONS_WLSAPPROXIMANT_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_WLSAPPROXIMANT_HPP_

/**
 * @file
 * Implementation of the function, obtained by WLS approximation.
 */

#include "WLSApproximant_fwd.hpp"
#include "medusa/bits/utils/numutils.hpp"

namespace mm {

template <typename basis_t>
WLSApproximant<basis_t>::WLSApproximant(const basis_t& basis, const vector_t& point,
               const std::vector<vector_t>& support, scalar_t scale,
               const Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>& coefficients,
               scalar_t residual) :
        basis_(basis), point_(point), support_(support), scale_(scale),
        coefficients_(coefficients), residual_(residual) {
    int n = support_.size();
    for (int i = 0; i < n; ++i) {
        support_[i] -= point;
        support_[i] /= scale;
    }
}

template <typename basis_t>
typename WLSApproximant<basis_t>::scalar_t
WLSApproximant<basis_t>::operator()(const vector_t& point) const {
    int m = basis_.size();
    vector_t local_point = (point - point_) / scale_;
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> b(m);
    for (int i = 0; i < m; ++i) {
        b(i) = basis_.eval(i, local_point, support_);
    }
    return b.dot(coefficients_);
}

template <typename basis_t>
template <typename operator_t>
typename WLSApproximant<basis_t>::scalar_t
WLSApproximant<basis_t>::operator()(const vector_t& point, const operator_t& op) const {
    int m = basis_.size();
    vector_t local_point = (point - point_) / scale_;
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> b(m);
    for (int i = 0; i < m; ++i) {
        b(i) = basis_.evalOp(i, local_point, op, support_, scale_);
    }
    return b.dot(coefficients_);
}


}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_WLSAPPROXIMANT_HPP_
