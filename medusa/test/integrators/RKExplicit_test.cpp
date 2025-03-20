#include "medusa/bits/integrators/RKExplicit.hpp"

#include "gtest/gtest.h"

namespace mm {

TEST(Integrators, RKExplicitStiff) {
    /*
     * Solving:
     * y' = -20(y-2), y(0) = 3
     *
     * with solution x = exp(-t)
     */
    std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)> func =
            [](double, const Eigen::VectorXd& y) {
                return -20*(y-2*Eigen::VectorXd::Ones(1));
            };

    Eigen::VectorXd y0(1); y0 << 3.0;
    double tmax = 10.0;
    double step = 0.05;
    auto solver_rk4 = integrators::Explicit::RK4().solve(func, 0.0, tmax, step, y0);
    auto solver_rk38 = integrators::Explicit::RK38().solve(func, 0.0, tmax, step, y0);
    auto solver_rk3 = integrators::Explicit::RK3().solve(func, 0.0, tmax, step, y0);
    auto solver_midpoint = integrators::Explicit::Midpoint().solve(func, 0.0, tmax, step, y0);
    auto solver_euler = integrators::Explicit::Euler().solve(func, 0.0, tmax, step, y0);

    auto stepper_rk4 = solver_rk4.begin();
    auto stepper_rk38 = solver_rk38.begin();
    auto stepper_rk3 = solver_rk3.begin();
    auto stepper_midpoint = solver_midpoint.begin();
    auto stepper_euler = solver_euler.begin();

    while (stepper_rk4) {
        ++stepper_rk4;
        ++stepper_rk38;
        ++stepper_rk3;
        ++stepper_midpoint;
        ++stepper_euler;
    }

    double correct = 2.0 + std::exp(-20.0*stepper_rk4.time());
    EXPECT_NEAR(correct, stepper_rk4.value()[0], 1e-9);
    EXPECT_NEAR(correct, stepper_rk38.value()[0], 1e-9);
    EXPECT_NEAR(correct, stepper_rk3.value()[0], 1e-7);
    EXPECT_NEAR(correct, stepper_midpoint.value()[0], 1e-6);
    EXPECT_NEAR(correct, stepper_euler.value()[0], 1e-4);
}

TEST(Integrators, RKExplicitCircle) {
    /*
     * Solving:
     * x' = -y
     * y' = x
     *
     * with solution (x, y) = (cos(t), sin(t))
     */

    std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)> func =
            [](double, const Eigen::VectorXd& y) {
                Eigen::VectorXd r(2);
                r(0) = -y(1);
                r(1) = y(0);
                return r;
            };

    double tmax = 2*PI;
    double step = 0.1;
    Eigen::VectorXd y0(2); y0 << 1.0, 0.0;
    auto solver_rk4 = integrators::Explicit::RK4().solve(func, 0.0, tmax, step, y0);
    auto solver_rk38 = integrators::Explicit::RK38().solve(func, 0.0, tmax, step, y0);
    auto solver_rk3 = integrators::Explicit::RK3().solve(func, 0.0, tmax, step, y0);
    auto solver_midpoint = integrators::Explicit::Midpoint().solve(func, 0.0, tmax, step, y0);
    auto solver_euler = integrators::Explicit::Euler().solve(func, 0.0, tmax, step, y0);

    auto stepper_rk4 = solver_rk4.begin();
    auto stepper_rk38 = solver_rk38.begin();
    auto stepper_rk3 = solver_rk3.begin();
    auto stepper_midpoint = solver_midpoint.begin();
    auto stepper_euler = solver_euler.begin();

    while (stepper_rk4) {
        ++stepper_rk4;
        ++stepper_rk38;
        ++stepper_rk3;
        ++stepper_midpoint;
        ++stepper_euler;
    }

    double correctx = std::cos(stepper_rk4.time()), correcty = std::sin(stepper_rk4.time());
    EXPECT_NEAR(correctx, stepper_rk4.value()[0], 1e-5);
    EXPECT_NEAR(correcty, stepper_rk4.value()[1], 1e-5);
    EXPECT_NEAR(correctx, stepper_rk38.value()[0], 1e-5);
    EXPECT_NEAR(correcty, stepper_rk38.value()[1], 1e-5);
    EXPECT_NEAR(correctx, stepper_rk3.value()[0], 1e-3);
    EXPECT_NEAR(correcty, stepper_rk3.value()[1], 1e-3);
    EXPECT_NEAR(correctx, stepper_midpoint.value()[0], 2e-2);
    EXPECT_NEAR(correcty, stepper_midpoint.value()[1], 2e-2);
    EXPECT_NEAR(correctx, stepper_euler.value()[0], 1e-0);
    EXPECT_NEAR(correcty, stepper_euler.value()[1], 1e-0);
}

TEST(Integrators, ExplicitOfOrderRK) {
    // Test that this compiles
    integrators::Explicit::of_order<1>();
    integrators::Explicit::of_order<2>();
    integrators::Explicit::of_order<3>();
    integrators::Explicit::of_order<4>();
}

TEST(Integrators, DISABLED_RungeKuttaUsageExample) {
    /// [Runge Kutta usage example]
    std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)> func =
            [](double, const Eigen::VectorXd& y) {
                return -y;
            };
    Eigen::VectorXd y0(1);
    y0 << 1.0;
    double tmax = 10.0;
    double dt = 0.1;
    auto integrator = integrators::Explicit::RK4().solve(func, 0.0, tmax, dt, y0);
    std::cout << integrator << std::endl;
    for (auto& step : integrator) {
        // You can use read write access to:
        step.value();
        step.time();
        // and check if this is the last step:
        step.is_last();
    }

    // Aditionally, one can iterate manually
    auto step = integrator.begin();
    while (step) {
        // do something with step
        ++step;
    }
    Eigen::VectorXd value = step.value();  // do something with value
    /// [Runge Kutta usage example]
}

}  // namespace mm
