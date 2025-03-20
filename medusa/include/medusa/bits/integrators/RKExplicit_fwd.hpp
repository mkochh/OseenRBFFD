#ifndef MEDUSA_BITS_INTEGRATORS_RKEXPLICIT_FWD_HPP_
#define MEDUSA_BITS_INTEGRATORS_RKEXPLICIT_FWD_HPP_

/**
 * @file
 * Runge-Kutta integrator declarations.
 *
 * @example test/integrators/RKExplicit_test.cpp
 */

#include <medusa/Config.hpp>
#include <Eigen/Core>
#include <medusa/bits/utils/numutils.hpp>
#include "RKExplicit_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>

namespace mm {

template <typename Scalar, int num_steps>
class AdamsBashforth;

/**
 * Class representing an explicit Runge-Kutta method.
 * @tparam Scalar type of scalars to use for calculations, e.g.\ double
 * @tparam num_stages number of stages of the method
 *
 * Usage example:
 * @snippet integrators/RKExplicit_test.cpp Runge Kutta usage example
 * @ingroup integrators
 */
template <typename Scalar, int num_stages>
class RKExplicit {
  public:
    typedef Scalar scalar_t;  ///< Floating point type used in the computations.
    enum { /** number of stages */ stages = num_stages };

  private:
    Eigen::Matrix<scalar_t, stages, 1> alpha_;  ///< Left row of the tableau.
    Eigen::Matrix<scalar_t, stages, stages> beta_;  ///< Central matrix of the tableau.
    Eigen::Matrix<scalar_t, stages, 1> gamma_;  ///< Bottom row of the tableau.

  public:
    /**
     * Construct a Runge-Kutta method with tableau
     * @f$
     * \begin{array}{l|l}
     * \alpha & \beta \\ \hline
     *        & \gamma^\mathsf{T}
     * \end{array}
     * @f$
     * @param alpha Vector of size `stages`.
     * @param beta Matrix of size `stages` times `stages`. Only strictly lower triangular part is
     * used.
     * @param gamma Vector of size `stages`.
     */
    RKExplicit(const Eigen::Matrix<scalar_t, stages, 1>& alpha,
               const Eigen::Matrix<scalar_t, stages, stages>& beta,
               const Eigen::Matrix<scalar_t, stages, 1>& gamma)
            : alpha_(alpha), beta_(beta), gamma_(gamma) {}

    /**
     * Integrator class for RK methods.
     */
    template <typename func_t>
    class Integrator {
        const RKExplicit<scalar_t, num_stages> method_;  ///< Method used in this integrator.
        scalar_t t0_;  ///< Start time.
        scalar_t dt_;  ///< Time step.
        int max_steps;  ///< Number of steps.
        Eigen::VectorXd y0_;  ///< Initial value.
        const func_t& func_;  ///< Function to integrate.

      public:
        /**
         * Constructs stepper. It is recommended to use the RKExplicit::solve() factory method
         * to construct this object.
         */
        Integrator(const RKExplicit<scalar_t, num_stages>& method, const func_t& func, scalar_t t0,
                   scalar_t t_max, scalar_t dt, Eigen::VectorXd y0) :
                method_(method), t0_(t0), dt_(dt), max_steps(iceil((t_max - t0) / dt)),
                y0_(std::move(y0)), func_(func) {}

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
            Eigen::VectorXd y;  ///< current value
            int cur_step;  ///< current step

          public:
            /// Construct an iterator at the initial values of the equation.
            explicit IterationStep(const Integrator& integrator)
                    : integrator(integrator), t(integrator.t0_), y(integrator.y0_), cur_step(0) {}

            /// Construct an (invalid) iterator pointing past the last step
            IterationStep(const Integrator& integrator, int) : integrator(integrator), t(), y(),
                                                               cur_step(integrator.max_steps + 1) {}

            /// Advance the stepper for one time step, returning the new value.
            IterationStep& operator++() {
                y = integrator.method_.step(integrator.func_, t, y, integrator.dt_);
                t += integrator.dt_;
                cur_step++;
                return *this;
            }

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

            /// Noop, used to conform to `std::iterator` requirements.
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
            Eigen::VectorXd value() const { return y; }

            /// Read-write access to current time.
            scalar_t& time() { return t; }

            /// Read-write access to current value.
            Eigen::VectorXd& value() { return y; }

            /// Output current state of the stepper.
            friend std::ostream& operator<<(std::ostream& os, const IterationStep& step) {
                return os << "t = " << step.t << "; y = " << step.y;
            }
        };

        typedef IterationStep iterator;  ///< iterator type
        typedef IterationStep const_iterator;  ///< const iterator type

        /// Creates stepper at positioned at initial value.
        iterator begin() { return IterationStep(*this); }

        /// Same as Integrator::begin()
        const_iterator cbegin() { return IterationStep(*this); }

        /// Creates an invalid stepper at positioned past the last step, meant to be used for
        /// comparison `it == end()` only
        iterator end() { return IterationStep(*this, 0); }

        /// Same as Integrator::end()
        const_iterator cend() { return IterationStep(*this, 0); }

        /// Output information about this integrator.
        friend std::ostream& operator<<(std::ostream& os, const Integrator& integrator) {
            return os << "Explicit RK integrator using " << integrator.method_ << " method from "
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
        return Integrator<func_t>(*this, func, t0, t_max, dt, y0);
    }

  private:
    /// Returns next value.
    template <typename func_t>
    Eigen::VectorXd step(const func_t& func, scalar_t t, const Eigen::VectorXd& y,
                         scalar_t dt) const {
        Eigen::Matrix<scalar_t, Eigen::Dynamic, stages> k(y.size(), static_cast<int>(stages));
        k.col(0) = func(t + alpha_(0) * dt, y);
        for (int i = 1; i < stages; ++i) {
            k.col(i) = func(t + alpha_(i) * dt,
                            y + dt * k.leftCols(i) * beta_.row(i).head(i).transpose());
        }
        return y + dt * k * gamma_;
    }

  public:
    /// Output the method's tableau for debugging.
    template <typename scalar_t, int stages>
    friend std::ostream& operator<<(std::ostream& os, const RKExplicit<scalar_t, stages>&);

    /// Implements Adams Bashforth multistep methods.
    template <typename Scalar2, int num_steps2>
    friend class AdamsBashforth;
};

/// Output the method's tableau for debugging.
template <typename scalar_t, int num_stages>
std::ostream& operator<<(std::ostream& os, const RKExplicit<scalar_t, num_stages>& method) {
    os << "RKExplicit with " << num_stages << " stages and "
       << "alpha = " << method.alpha_ << ", "
       << "beta = " << method.beta_ << ", "
       << "gamma = " << method.gamma_;
    return os;
}

/**
 * Namespace containing most known methods for integrating ODEs.
 */
namespace integrators {

/**
 * Namespace containing factory functions for explicit single step integrators.
 * Factory functions are provided for the most common ones, such as Euler's method, Midpoint rule,
 * RK4, ...
 *
 * Usage example:
 * @snippet RKExplicit_test.cpp Runge Kutta usage example
 */
namespace Explicit {

/// Standard RK4 4th order method
template <class scalar_t = double>
RKExplicit<scalar_t, 4> RK4();

/// 3/8 rule 4th order method
template <class scalar_t = double>
RKExplicit<scalar_t, 4> RK38();

/// Explicit Euler's method, 1st order
template <class scalar_t = double>
RKExplicit<scalar_t, 1> Euler();

/// Explicit midpoint rule, 2nd order
template <class scalar_t = double>
RKExplicit<scalar_t, 2> Midpoint();

/// Kutta's third order method
template <class scalar_t = double>
RKExplicit<scalar_t, 3> RK3();

/// Fifth order method appearing in Fehlberg's method
template <class scalar_t = double>
RKExplicit<scalar_t, 6> Fehlberg5();

/// @cond
namespace of_order_internal {
template <int order, class scalar_t>
struct _of_order_impl { static RKExplicit<scalar_t, order> of_order(); };
}
/// @endcond

/**
 * Returns Runge Kutta explicit method of requested order with given floating point type.
 */
template <int order, class scalar_t = double>
static RKExplicit<scalar_t, order> of_order() {
    return of_order_internal::_of_order_impl<order, scalar_t>::of_order();
}

}  // namespace Explicit
}  // namespace integrators
}  // namespace mm

#endif  // MEDUSA_BITS_INTEGRATORS_RKEXPLICIT_FWD_HPP_
