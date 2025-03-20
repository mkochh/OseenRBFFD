#ifndef MEDUSA_BITS_INTEGRATORS_ADAMSBASHFORTH_HPP_
#define MEDUSA_BITS_INTEGRATORS_ADAMSBASHFORTH_HPP_

/**
 * @file
 * Adams-Bashforth integrator implementation.
 */

#include "AdamsBashforth_fwd.hpp"
#include "RKExplicit.hpp"

namespace mm {

template <typename Scalar, int num_steps>
template <typename func_t, typename initial_method_t>
typename AdamsBashforth<Scalar, num_steps>::template
        Integrator<func_t, initial_method_t>::IterationStep&
AdamsBashforth<Scalar, num_steps>::Integrator<func_t,
                                              initial_method_t>::IterationStep::operator++() {
    cur_step++;
    Eigen::VectorXd next;
    if (cur_step < steps) {  // don't have enough steps, do a RK method.
        next = integrator.initial_method_.step(
                integrator.func_, t, last_ys.col(steps - 1), integrator.dt_);
    } else {
        next = integrator.method_.step(integrator.func_, t, last_ys, integrator.dt_);
    }
    t += integrator.dt_;
    // shift values
    for (int i = 0; i < steps-1; ++i) {
        last_ys.col(i) = last_ys.col(i+1);
    }
    last_ys.col(steps-1) = next;
    return *this;
}

/// Output the method's tableau for debugging.
template <typename scalar_t, int num_steps>
std::ostream& operator<<(std::ostream& os, const AdamsBashforth<scalar_t, num_steps>& method) {
    os << "AdamsBashforth with " << num_steps << " stages and "
       << "b = " << method.b_.transpose();
    return os;
}

namespace integrators {
namespace ExplicitMultistep {
/// @cond
namespace of_order_internal {

template <int order, class scalar_t>
AdamsBashforth<scalar_t, order> _of_order_impl<order, scalar_t>::of_order() {
    static_assert(order > 0, "Order must be positive.");
    static_assert(order < 6, "Methods of so high orders are not implemented yet, but it's"
                             " very simple to add them!");
    throw;
}

/// Specialization of order 1
template <class scalar_t>
struct _of_order_impl<1, scalar_t> {
    static AdamsBashforth<scalar_t, 1> of_order() { return AB1(); }
};

/// Specialization of order 2
template <class scalar_t>
struct _of_order_impl<2, scalar_t> {
    static AdamsBashforth<scalar_t, 2> of_order() { return AB2(); }
};

/// Specialization of order 3
template <class scalar_t>
struct _of_order_impl<3, scalar_t> {
    static AdamsBashforth<scalar_t, 3> of_order() { return AB3(); }
};

/// Specialization of order 4
template <class scalar_t>
struct _of_order_impl<4, scalar_t> {
    static AdamsBashforth<scalar_t, 4> of_order() { return AB4(); }
};

/// Specialization of order 5
template <class scalar_t>
struct _of_order_impl<5, scalar_t> {
    static AdamsBashforth<scalar_t, 5> of_order() { return AB5(); }
};

}  // namespace of_order_internal
/// @endcond
}  // namespace ExplicitMultistep
}  // namespace integrators
}  // namespace mm

#endif  // MEDUSA_BITS_INTEGRATORS_ADAMSBASHFORTH_HPP_
