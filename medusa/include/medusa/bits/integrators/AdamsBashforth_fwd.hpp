#ifndef MEDUSA_BITS_INTEGRATORS_ADAMSBASHFORTH_FWD_HPP_
#define MEDUSA_BITS_INTEGRATORS_ADAMSBASHFORTH_FWD_HPP_

/**
 * @file
 * Adams-Bashforth integrator declarations.
 *
 * @example test/integrators/AdamsBashforth_test.cpp
 */

#include <medusa/Config.hpp>
#include <Eigen/Core>
#include <medusa/bits/utils/numutils.hpp>
#include "RKExplicit_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>

namespace mm {

/**
 * Class representing an AdamsBashforth method, an explicit linear multistep method.
 * @tparam Scalar Type of scalars to use for calculations, e.g.\ `double`.
 * @tparam num_steps Number of previous steps the method uses.
 *
 * Usage example:
 * @snippet integrators/AdamsBashforth_test.cpp Adams Bashforth usage example
 * @ingroup integrators
 */
template <typename Scalar, int num_steps>
class AdamsBashforth {
  public:
    typedef Scalar scalar_t;  ///< floating point type
    enum {
        steps = num_steps  ///< number of steps of the multistep method
    };
  private:
    Eigen::Matrix<scalar_t, steps, 1> b_;  ///< coefficients of the multistep method

  public:
    /**
     * Construct an `num_steps`-step Adams Bashforth method with coefficients `b`, given as
     * @f$y_{n+s} = y_{n+s-1} + h \sum_{j=0}^{s-1} b_j f(t_j, y_j)@f$.
     * @param b Vector of size `num_steps` representing the coefficients of the method.
     */
    AdamsBashforth(const Eigen::Matrix<scalar_t, steps, 1>& b) : b_(b) {
        assert_msg(std::abs(b_.sum() - 1.0) < 1e-15, "AB method coefficients do not sum up to 1, "
                                                     "but %.16f instead.", b_.sum());
    }

    /**
     * Integrator class for AB methods.
     * @tparam starting_method_t Types of the underlying starting method, required to get first
     * `num_steps` values.
     * @warning
     * You should choose a method that provides same order of accuracy, is a single step method and
     * has a fixed time step. A default choice is the `num_steps` stage Runge Kutta method,
     * which might not always be correct. Predefined methods that are called by name, such
     * as ExplicitMultistep::AB5(), have correct initial methods.
     */
    template <typename func_t, typename initial_method_t = RKExplicit<Scalar, num_steps>>
    class Integrator {
        const initial_method_t initial_method_;  ///< method used for first `num_steps` steps
        const AdamsBashforth<scalar_t, num_steps> method_;  ///< method used for all next steps
        scalar_t t0_;  ///< start time
        scalar_t dt_;  ///< time step
        int max_steps;  ///< number of steps
        Eigen::VectorXd y0_;  ///< initial value
        const func_t& func_;  ///< function to integrate

      public:
        /**
         * Constructs stepper with given initial method. It is recommended to use the
         * AdamsBashforth::solve() factory method to construct this object.
         */
        Integrator(const AdamsBashforth<scalar_t, num_steps>& method,
                   const initial_method_t& initial_method,
                   const func_t& func, scalar_t t0, scalar_t t_max, scalar_t dt, Eigen::VectorXd y0)
                : initial_method_(initial_method), method_(method), t0_(t0),
                  dt_(dt), max_steps(iceil((t_max - t0) / dt)), y0_(std::move(y0)),
                  func_(func) {}

        /**
         * Class representing a step in the integration process. This class satisfies
         * `std::forward_iterator` requirements and can therefore be used with STL-type algorithms.
         */
        class IterationStep : public std::iterator<
                std::forward_iterator_tag,   // iterator_category
                IterationStep,               // value_type
                int,                         // difference_type
                const IterationStep*,        // pointer
                IterationStep                // reference
        > {
            const Integrator& integrator;  ///< reference to underlying integrator
            scalar_t t;  ///< current time
            Eigen::Matrix<scalar_t, Eigen::Dynamic, steps> last_ys;  ///< current value
            int cur_step;  ///< current step

          public:
            /// Construct an iterator at the initial values of the equation.
            explicit IterationStep(const Integrator& integrator)
                    : integrator(integrator), t(integrator.t0_),
                      last_ys(integrator.y0_.size(), static_cast<int>(steps)), cur_step(0) {
                last_ys.col(steps-1) = integrator.y0_;
            }

            /// Construct an (invalid) iterator pointing past the last step
            IterationStep(const Integrator& integrator, int) :
                    integrator(integrator), t(), last_ys(), cur_step(integrator.max_steps + 1) {}

            /// Advance the stepper for one time step, returning the new value.
            IterationStep& operator++();

            /**
             * Advance the stepper for one time step, returning the old value.
             * @warning usually the prefix increment `++step` is preferred, due to not having a
             * temporary variable.
             */
            IterationStep operator++(int) {
                IterationStep retval = *this;
                ++(*this);
                return retval;
            }

            /// Compare two steppers if they are on the same step.
            bool operator==(const IterationStep& other) const { return cur_step == other.cur_step; }

            /// Negation of IterationStep::operator==.
            bool operator!=(const IterationStep& other) const { return !(*this == other); }

            /// Noop, used to conform to `std::iterator` requirements
            IterationStep& operator*() { return *this; }

            /// const version of IterationStep::operator*
            const IterationStep& operator*() const { return *this; }

            /// Returns `false` if integrator went past the last step and `true` otherwise.
            explicit operator bool() {
                return cur_step <= integrator.max_steps;
            }

            /// Returns true if integrator just completed its last step.
            bool is_last() const {
                return cur_step == integrator.max_steps;
            }

            /// Get current time.
            scalar_t time() const { return t; }

            /// Get current value.
            typename Eigen::Matrix<scalar_t, Eigen::Dynamic, steps>::ConstColXpr value() const {
                return last_ys.col(steps-1);
            }

            /// Read-write access to current time
            scalar_t& time() { return t; }

            /// Read-write access to current value.
            typename Eigen::Matrix<scalar_t, Eigen::Dynamic, steps>::ColXpr value() {
                return last_ys.col(steps-1);
            }

            /// Output current state of the stepper
            friend std::ostream& operator<<(std::ostream& os, const IterationStep& step) {
                return os << "t = " << step.t << "; y = " << step.last_ys.col(steps-1);
            }
        };

        typedef IterationStep iterator;  ///< iterator type
        typedef IterationStep const_iterator;  ///< const iterator type

        /// Creates stepper at positioned at initial value.
        iterator begin() { return IterationStep(*this); }

        /// Same as Integrator::begin().
        const_iterator cbegin() { return IterationStep(*this); }

        /// Creates an invalid stepper at positioned past the last step, meant to be used for
        /// comparison `it == end()` only.
        iterator end() { return IterationStep(*this, 0); }

        /// Same as Integrator::end().
        const_iterator cend() { return IterationStep(*this, 0); }

        /// Output information about this integrator.
        friend std::ostream& operator<<(std::ostream& os, const Integrator& integrator) {
            return os << "Multistep integrator using " << integrator.method_ << " method from "
                      << "time " << integrator.t0_ << " with step " << integrator.dt_ << " for "
                      << integrator.max_steps << " steps.";
        }
    };

    /**
     * Returns a solver using this method.
     * @param func Function to integrate.
     * @param t0 Start time.
     * @param t_max End time.
     * @param dt Time step.
     * @param y0 Initial value.
     * @warning Last step might not finish exactly on `t_max`. The number of steps is calculated as
     * `std::ceil((t_max - t0)/dt)`.
     * @return A solver object.
     */
    template <typename func_t>
    Integrator<func_t> solve(const func_t& func, scalar_t t0, scalar_t t_max, scalar_t dt,
                             const Eigen::VectorXd& y0) const {
        return Integrator<func_t, RKExplicit<scalar_t, num_steps>>(
                *this, integrators::Explicit::of_order<num_steps>(), func, t0, t_max, dt, y0);
    }

    /// Same as AdamsBashforth::solve, but offers the option of supplying your own initial method.
    template <typename initial_method_t, typename func_t>
    Integrator<func_t, initial_method_t> solve(const initial_method_t& method, const func_t& func,
                                               scalar_t t0, scalar_t t_max, scalar_t dt,
                                               const Eigen::VectorXd& y0) const {
        return Integrator<func_t, initial_method_t>(*this, method, func, t0, t_max, dt, y0);
    }

  private:
    /// Returns next value
    template <typename func_t>
    Eigen::VectorXd step(const func_t& func, scalar_t t,
                         const Eigen::Matrix<scalar_t, Eigen::Dynamic, steps>& last_ys,
                         scalar_t dt) const {
        Eigen::VectorXd result = last_ys.col(steps-1);
        for (int i = 0; i < steps; ++i) {
            result += dt * b_(i) * func(t - i*dt, last_ys.col(i));
        }
        return result;
    }

  public:
    /// Output the method's tableau for debugging.
    template <typename scalar_t, int stages>
    friend std::ostream& operator<<(std::ostream& os, const AdamsBashforth<scalar_t, stages>&);
};

namespace integrators {

/**
 * Namespace containing factory functions for explicit linear multistep integrators.
 *
 * Usage example:
 * @snippet AdamsBashforth_test.cpp Adams Bashforth usage example
 */
namespace ExplicitMultistep {

/// Standard Euler's method.
template <class scalar_t = double>
static AdamsBashforth<scalar_t, 1> AB1() {
    Eigen::Matrix<scalar_t, 1, 1> b(1); b << 1.0;
    return {b};
}

/// Two step AB method.
template <class scalar_t = double>
static AdamsBashforth<scalar_t, 2> AB2() {
    Eigen::Matrix<scalar_t, 2, 1> b(2); b << -1./2., 3./2.;
    return {b};
}

/// Three step AB method.
template <class scalar_t = double>
static AdamsBashforth<scalar_t, 3> AB3() {
    Eigen::Matrix<scalar_t, 3, 1> b(3); b << 5./12., -4./3., 23./12.;
    return {b};
}

/// Four step AB method.
template <class scalar_t = double>
static AdamsBashforth<scalar_t, 4> AB4() {
    Eigen::Matrix<scalar_t, 4, 1> b(4); b << -3./8., 37./24., -59./24., 55./24.;
    return {b};
}

/// Five step AB method.
template <class scalar_t = double>
static AdamsBashforth<scalar_t, 5> AB5() {
    Eigen::Matrix<scalar_t, 5, 1> b(5);
    b << 251./720., -637./360., 109./30., -1387./360., 1901./720.;
    return {b};
}

// Internal namespace containing implementation details for Explicit::of_order function.
/// @cond
namespace of_order_internal {
template <int order, class scalar_t>
struct _of_order_impl {
    static AdamsBashforth<scalar_t, order> of_order();
};
}  // namespace of_order_internal
/// @endcond

/**
 * Returns Adams-Bashforth explicit method of requested order with given floating point type.
 */
template <int order, class scalar_t = double>
static AdamsBashforth<scalar_t, order> of_order() {
    return of_order_internal::_of_order_impl<order, scalar_t>::of_order();
}

}  // namespace ExplicitMultistep
}  // namespace integrators
}  // namespace mm

#endif  // MEDUSA_BITS_INTEGRATORS_ADAMSBASHFORTH_FWD_HPP_
