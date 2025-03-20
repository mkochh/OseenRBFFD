#ifndef MEDUSA_BITS_APPROXIMATIONS_WLS_HPP_
#define MEDUSA_BITS_APPROXIMATIONS_WLS_HPP_

/**
 * @file
 * Implementation of weighted least squares approximation.
 */

#include "WLS_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/types/Range.hpp>
#include <medusa/bits/utils/numutils.hpp>
#include "JacobiSVDWrapper.hpp"
#include "WLSApproximant.hpp"

namespace mm {

template <class basis_t, class weight_t, class scale_t, class solver_t>
void WLS<basis_t, weight_t, scale_t, solver_t>::compute(const vector_t& point,
                                                        const std::vector<vector_t>& support) {
    int n = support.size();
    assert_msg(basis_.size() <= n,
               "WLS cannot have more basis functions than points in support domain, but %d "
               "support points and %d basis functions were given.", n, basis_.size());
    assert_msg(basis_.size() > 0, "Cannot construct an approximation with an empty basis.");
    scale_ = scale_t::scale(point, support);
    local_coordinates_.resize(n);
    for (int i = 0; i < n; ++i) {
        local_coordinates_[i] = (support[i] - point) / scale_;
    }
    point_ = point;
    solver_.compute(getMatrix(local_coordinates_));
}

template <class basis_t, class weight_t, class scale_t, class solver_t>
typename WLS<basis_t, weight_t, scale_t, solver_t>::ei_matrix_t
WLS<basis_t, weight_t, scale_t, solver_t>::getMatrix(
        const std::vector<vector_t>& local_coordinates) const {
    // Evaluate basis functions in given points
    int n = local_coordinates.size();
    int m = basis_.size();
    ei_matrix_t WB(m, n);
    for (int j = 0; j < n; j++) {
        scalar_t w = weight_(local_coordinates[j]);
        for (int i = 0; i < m; i++) {
            WB(i, j) = w * basis_.eval(i, local_coordinates[j], local_coordinates);
        }
    }
    return WB;
}

template <class basis_t, class weight_t, class scale_t, class solver_t>
Eigen::Matrix<typename basis_t::scalar_t, Eigen::Dynamic, 1>
WLS<basis_t, weight_t, scale_t, solver_t>::getShape() const {
    int m = basis_.size();
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> b(m);
    for (int i = 0; i < m; i++) {
        b(i) = basis_.evalAt0(i, local_coordinates_);
    }
    return getShapeFromRhs(b);
}

/// @cond
template <class basis_t, class weight_t, class scale_t, class solver_t>
template <class operator_t>
Eigen::Matrix<typename basis_t::scalar_t, Eigen::Dynamic, 1>
WLS<basis_t, weight_t, scale_t, solver_t>::getShape(const operator_t& op) const {
    int m = basis_.size();
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> Lb(m);
    for (int i = 0; i < m; i++) {
        Lb(i) = basis_.evalOpAt0(i, op, local_coordinates_, scale_);
    }
    return getShapeFromRhs(Lb);
}

template <class basis_t, class weight_t, class scale_t, class solver_t>
template <class Derived>
Eigen::Matrix<typename basis_t::scalar_t, Eigen::Dynamic, 1>
WLS<basis_t, weight_t, scale_t, solver_t>::getShapeFromRhs(
        const Eigen::MatrixBase<Derived>& Lb) const {
    int n = local_coordinates_.size();
    // Solve (WB)^T x = Lb
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> result = solver_.solve(Lb);
    // Multiply with weight
    for (int j = 0; j < n; j++) {
        result(j) *= weight_(local_coordinates_[j]);
    }
    return result;
}

template <class basis_t, class weight_t, class scale_t, class solver_t>
Eigen::Matrix<typename basis_t::scalar_t, Eigen::Dynamic, 1>
WLS<basis_t, weight_t, scale_t, solver_t>::getWeights() const {
    int n = local_coordinates_.size();
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> w(n);
    for (int j = 0; j < n; j++) {
        w(j) = weight_(local_coordinates_[j]);
    }
    return w;
}
/// @endcond

/// Output basic info about given WLS class.
template <class basis_t, class weight_t, class scale_t, class solver_t>
std::ostream& operator<<(std::ostream& os, const WLS<basis_t, weight_t, scale_t, solver_t>& wls) {
    os << "WLS:\n"
       << "    dimension: " << wls.dim << '\n'
       << "    basis size: " << wls.basis_.size() << '\n'
       << "    basis: "  << wls.basis_ << '\n'
       << "    weight: "  << wls.weight_ << '\n'
       << "    scale: "  << scale_t() << '\n';
    if (wls.point_[0] != wls.point_[0]) {
        os << "    last point: not used yet";
    } else {
        os << "    last point: " << wls.point_;
    }
    return os;
}

template <class basis_t, class weight_t, class scale_t, class solver_t>
template <typename Derived>
WLSApproximant<basis_t>
WLS<basis_t, weight_t, scale_t, solver_t>::getApproximant(
        const vector_t& point, const std::vector<vector_t>& support,
        const Eigen::MatrixBase<Derived>& values) const {
    int n = support.size();
    assert_msg(values.size() == n, "Number of given values %d must be equal to support size %d.",
               values.size(), n);
    assert_msg(basis_.size() <= n,
               "WLS cannot have more basis functions than points in support domain, but %d "
               "support points and %d basis functions were given.", n, basis_.size());
    double scale = scale_t::scale(point, support);
    std::vector<vector_t> local_coordinates(n);
    for (int i = 0; i < n; ++i) {
        local_coordinates[i] = (support[i] - point) / scale;
    }
    solver_t solver;
    ei_matrix_t matrix = getMatrix(local_coordinates).transpose();
    solver.compute(matrix);

    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> w(n);
    for (int j = 0; j < n; j++) {
        w(j) = weight_(local_coordinates[j]);
    }
    Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> coefficients = solver.solve(values.cwiseProduct(w));
    scalar_t residual = (matrix*coefficients - values).squaredNorm();
    return WLSApproximant<basis_t>(basis_, point, support, scale, coefficients, residual);
}

}  // namespace mm

#endif  // MEDUSA_BITS_APPROXIMATIONS_WLS_HPP_
