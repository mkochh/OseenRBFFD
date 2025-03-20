#ifndef MEDUSA_BITS_INTEGRATORS_RKEXPLICIT_HPP_
#define MEDUSA_BITS_INTEGRATORS_RKEXPLICIT_HPP_

/**
 * @file
 * Runge-Kutta integrator implementation.
 */

#include "RKExplicit_fwd.hpp"

namespace mm {
namespace integrators {
namespace Explicit {

template <class scalar_t>
RKExplicit<scalar_t, 4> RK4() {
    Eigen::Matrix<scalar_t, 4, 1> alpha;
    alpha << 0, 0.5, 0.5, 1.0;
    Eigen::Matrix<scalar_t, 4, 4> beta;
    beta << 0.0, 0.0, 0.0, 0.0,
            0.5, 0.0, 0.0, 0.0,
            0.0, 0.5, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0;
    Eigen::Matrix<scalar_t, 4, 1> gamma;
    gamma << 1. / 6., 2. / 6., 2. / 6., 1. / 6.;
    return {alpha, beta, gamma};
}

template <class scalar_t>
RKExplicit<scalar_t, 4> RK38() {
    Eigen::Matrix<scalar_t, 4, 1> alpha;
    alpha << 0, 1. / 3., 2. / 3., 1.0;
    Eigen::Matrix<scalar_t, 4, 4> beta;
    beta << 0.0, 0.0, 0.0, 0.0,
            1. / 3., 0.0, 0.0, 0.0,
            -1. / 3., 1.0, 0.0, 0.0,
            1.0, -1.0, 1.0, 0.0;
    Eigen::Matrix<scalar_t, 4, 1> gamma;
    gamma << 1. / 8., 3. / 8., 3. / 8., 1. / 8.;
    return {alpha, beta, gamma};
}

template <class scalar_t>
RKExplicit<scalar_t, 1> Euler() {
    Eigen::Matrix<scalar_t, 1, 1> alpha;
    alpha << 0.0;
    Eigen::Matrix<scalar_t, 1, 1> beta;
    beta << 0.0;
    Eigen::Matrix<scalar_t, 1, 1> gamma;
    gamma << 1.0;
    return {alpha, beta, gamma};
}

template <class scalar_t>
RKExplicit<scalar_t, 2> Midpoint() {
    Eigen::Matrix<scalar_t, 2, 1> alpha;
    alpha << 0.0, 0.5;
    Eigen::Matrix<scalar_t, 2, 2> beta;
    beta << 0.0, 0.0, 0.5, 0.0;
    Eigen::Matrix<scalar_t, 2, 1> gamma;
    gamma << 0.0, 1.0;
    return {alpha, beta, gamma};
}

template <class scalar_t>
RKExplicit<scalar_t, 3> RK3() {
    Eigen::Matrix<scalar_t, 3, 1> alpha;
    alpha << 0, 0.5, 1.0;
    Eigen::Matrix<scalar_t, 3, 3> beta;
    beta << 0.0, 0.0, 0.0,
            0.5, 0.0, 0.0,
            -1.0, 2.0, 0.0;
    Eigen::Matrix<scalar_t, 3, 1> gamma;
    gamma << 1. / 6., 2. / 3., 1. / 6.;
    return {alpha, beta, gamma};
}

template <class scalar_t>
RKExplicit<scalar_t, 6> Fehlberg5() {
    Eigen::Matrix<scalar_t, 6, 1> alpha;
    alpha << 0, 1. / 4., 3. / 8., 12. / 13., 1., 1. / 2.;
    Eigen::Matrix<scalar_t, 6, 6> beta;
    beta << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            1. / 4., 0.0, 0.0, 0.0, 0.0, 0.0,
            3. / 32., 9. / 32., 0.0, 0.0, 0.0, 0.0,
            1932. / 2197, -7200. / 2197., 7296. / 2197., 0.0, 0.0, 0.0,
            439. / 216., -8., 3680. / 513., -845. / 4104., 0.0, 0.0,
            -8. / 27., 2., -3544. / 2565., 1859. / 4104., -11. / 40., 0.0;
    Eigen::Matrix<scalar_t, 6, 1> gamma;
    gamma << 16. / 135., 0.0, 6656. / 12825., 28561. / 56430., -9. / 50., 2. / 55.;
    return {alpha, beta, gamma};
}

/// @cond
namespace of_order_internal {

template <int order, class scalar_t>
RKExplicit<scalar_t, order> _of_order_impl<order, scalar_t>::of_order() {
    static_assert(order > 0, "Order must be positive.");
    static_assert(order < 5, "Methods of higher orders do not correspond with the number of "
                             "stages, you must call them manually. If you got this error "
                             "calling AB5() method, supply e.g. Fehlberg5 as the first "
                             "parameter and its type as the template parameter.");
    throw;  // to avoid warnings
}

/// Specialization of order 1, Euler's method.
template <class scalar_t>
struct _of_order_impl<1, scalar_t> {
    static RKExplicit<scalar_t, 1> of_order() { return Euler(); }
};

/// Specialization of order 2, Midpoint method.
template <class scalar_t>
struct _of_order_impl<2, scalar_t> {
    static RKExplicit<scalar_t, 2> of_order() { return Midpoint(); }
};

/// Specialization of order 3, RK3.
template <class scalar_t>
struct _of_order_impl<3, scalar_t> {
    static RKExplicit<scalar_t, 3> of_order() { return RK3(); }
};

/// Specialization of order 4, RK4.
template <class scalar_t>
struct _of_order_impl<4, scalar_t> {
    static RKExplicit<scalar_t, 4> of_order() { return RK4(); }
};

}  // namespace of_order_internal
/// @endcond

}  // namespace Explicit
}  // namespace integrators
}  // namespace mm

#endif  // MEDUSA_BITS_INTEGRATORS_RKEXPLICIT_HPP_
