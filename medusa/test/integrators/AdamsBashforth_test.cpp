#include "medusa/bits/integrators/AdamsBashforth.hpp"

#include "gtest/gtest.h"

namespace mm {

TEST(Integrators, ABStiff) {
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
    double step = 0.005;
    auto solver_ab1 = integrators::ExplicitMultistep::AB1().solve(func, 0.0, tmax, step, y0);
    auto solver_ab2 = integrators::ExplicitMultistep::AB2().solve(func, 0.0, tmax, step, y0);
    auto solver_ab3 = integrators::ExplicitMultistep::AB3().solve(func, 0.0, tmax, step, y0);
    auto solver_ab4 = integrators::ExplicitMultistep::AB4().solve(func, 0.0, tmax, step, y0);
    auto solver_ab5 = integrators::ExplicitMultistep::AB5().solve<RKExplicit<double, 6>>(
            integrators::Explicit::Fehlberg5(),
            func, 0.0, tmax, step, y0);

    auto stepper_ab1 = solver_ab1.begin();
    auto stepper_ab2 = solver_ab2.begin();
    auto stepper_ab3 = solver_ab3.begin();
    auto stepper_ab4 = solver_ab4.begin();
    auto stepper_ab5 = solver_ab5.begin();

    while (stepper_ab4) {
        ++stepper_ab1;
        ++stepper_ab2;
        ++stepper_ab3;
        ++stepper_ab4;
        ++stepper_ab5;
    }

    double correct = 2.0 + std::exp(-20.0*stepper_ab4.time());
    EXPECT_NEAR(correct, stepper_ab1.value()[0], 1e-4);
    EXPECT_NEAR(correct, stepper_ab2.value()[0], 1e-5);
    EXPECT_NEAR(correct, stepper_ab3.value()[0], 1e-6);
    EXPECT_NEAR(correct, stepper_ab4.value()[0], 1e-7);
    EXPECT_NEAR(correct, stepper_ab5.value()[0], 1e-8);
}

TEST(Integrators, ABCircle) {
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
    auto solver_ab1 = integrators::ExplicitMultistep::AB1().solve(func, 0.0, tmax, step, y0);
    auto solver_ab2 = integrators::ExplicitMultistep::AB2().solve(func, 0.0, tmax, step, y0);
    auto solver_ab3 = integrators::ExplicitMultistep::AB3().solve(func, 0.0, tmax, step, y0);
    auto solver_ab4 = integrators::ExplicitMultistep::AB4().solve(func, 0.0, tmax, step, y0);
    auto solver_ab5 = integrators::ExplicitMultistep::AB5().solve<RKExplicit<double, 6>>(
            integrators::Explicit::Fehlberg5(),
            func, 0.0, tmax, step, y0);

    auto stepper_ab1 = solver_ab1.begin();
    auto stepper_ab2 = solver_ab2.begin();
    auto stepper_ab3 = solver_ab3.begin();
    auto stepper_ab4 = solver_ab4.begin();
    auto stepper_ab5 = solver_ab5.begin();

    while (stepper_ab4) {
        ++stepper_ab1;
        ++stepper_ab2;
        ++stepper_ab3;
        ++stepper_ab4;
        ++stepper_ab5;
    }

    double correctx = std::cos(stepper_ab4.time()), correcty = std::sin(stepper_ab4.time());
    EXPECT_NEAR(correctx, stepper_ab4.value()[0], 5e-4);
    EXPECT_NEAR(correcty, stepper_ab4.value()[1], 5e-4);
    EXPECT_NEAR(correctx, stepper_ab3.value()[0], 5e-3);
    EXPECT_NEAR(correcty, stepper_ab3.value()[1], 5e-3);
    EXPECT_NEAR(correctx, stepper_ab2.value()[0], 5e-2);
    EXPECT_NEAR(correcty, stepper_ab2.value()[1], 5e-2);
    EXPECT_NEAR(correctx, stepper_ab1.value()[0], 1e-0);
    EXPECT_NEAR(correcty, stepper_ab1.value()[1], 1e-0);
}

TEST(Integrators, ExplicitOfOrderAB) {
    // Test that this compiles
    integrators::ExplicitMultistep::of_order<1>();
    integrators::ExplicitMultistep::of_order<2>();
    integrators::ExplicitMultistep::of_order<3>();
    integrators::ExplicitMultistep::of_order<4>();
    integrators::ExplicitMultistep::of_order<5>();
}

TEST(Integrators, DISABLED_ABUsageExample) {
    /// [Adams Bashforth usage example]
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
    double dt = 0.1;
    Eigen::VectorXd y0(2); y0 << 1.0, 0.0;
    auto solver = integrators::ExplicitMultistep::AB3().solve(func, 0.0, tmax, dt, y0);
    std::cout << solver << std::endl;
    for (auto step : solver) {
        std::cout << step << std::endl;
        if (step.is_last()) {
            std::cout << "This is the last step." << std::endl;
        }
        // Read/write access to current time and value.
        step.time();
        step.value();
    }
    /// [Adams Bashforth usage example]
}

}  // namespace mm
