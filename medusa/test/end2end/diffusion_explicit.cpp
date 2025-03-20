#include <medusa/Medusa_fwd.hpp>

#include "gtest/gtest.h"

namespace mm {

class DiffusionExplicitTest : public ::testing::Test {
  public:
    double dx;
    int n, m;
    double time, dt, t_steps;
    double sigma;
    DomainDiscretization<Vec2d> domain;
    Range<int> interior, boundary;
    RaggedShapeStorage<Vec2d, std::tuple<Lap<2>>> shapes;

    DiffusionExplicitTest() : dx(1. / 50.), n(12), m(2), time(0.05), dt(1e-5),
                              t_steps(std::ceil(time / dt)), sigma(1.0), domain(prep_domain()),
                              interior(domain.interior()), boundary(domain.boundary()),
                              shapes() {
        prep_shapes();
    }

  protected:
    DomainDiscretization<Vec2d> prep_domain() const {
        // prep domain
        BoxShape<Vec2d> b({0, 0}, {1, 1});
        auto d = b.discretizeWithStep(dx);
        d.findSupport(FindClosest(n));
        return d;
    }

    void prep_shapes() {
        // prep shape funcs
        WLS<Monomials<Vec2d>, GaussianWeight<Vec2d>, ScaleToClosest> wls(m, sigma);
        sh::operator_tuple<sh::lap, 2>::type operators;
        shapes.resize(domain.supportSizes());
        computeShapes(domain, wls, domain.all(), operators, &shapes);
    }

    static double analytical(const Vec2d& pos, double t, double a, double D, size_t N) {
        double T = 0;
        double f = PI / a;
        for (size_t n = 1; n < N; n = n + 2) {
            for (size_t m = 1; m < N; m = m + 2) {
                T += 16.0 / f / f / (n * m) * std::sin(n * f * pos[0]) * std::sin(m * f * pos[1]) *
                     std::exp(-D * t * ((n * n + m * m) * f * f));
            }
        }
        return T;
    }

    double LinfError(double t, const ScalarFieldd& value) {
        Range<double> E2(domain.size(), 0);
        for (auto& c : interior) {
            E2[c] = std::abs(value[c] - analytical(domain.pos(c), t, 1, 1, 50));
        }
        return *std::max_element(E2.begin(), E2.end());
    }

  public:
    double solveBasic() {
        ScalarFieldd T1(domain.size());
        T1[interior] = 1.0;
        T1[boundary] = 0.0;
        ScalarFieldd T2 = T1;
        int tt;
        for (tt = 0; tt < t_steps; ++tt) {
            // new temperature
            for (auto& c : interior) {
                double Lap = 0;
                for (int i = 0; i < n; ++i) Lap += shapes.laplace(c, i) * T1[domain.support(c, i)];
                T2[c] = T1[c] + dt * Lap;
            }
            // time advance
            T1.swap(T2);
        }
        return LinfError(tt * dt, T1);
    }

    double solveOperators() {
        ScalarFieldd T1(domain.size());
        T1[interior] = 1.0;
        T1[boundary] = 0.0;
        ScalarFieldd T2 = T1;

        auto op = shapes.explicitOperators();
        int tt;
        for (tt = 0; tt < t_steps; ++tt) {
            // new temp.
            for (auto& c : interior) {
                T2[c] = T1[c] + dt * op.lap(T1, c);
            }
            // time advance
            T1.swap(T2);
        }
        return LinfError(tt * dt, T1);
    }

    double solveRKEuler() {
        auto op = shapes.explicitOperators();

        auto dv_dt = [&](double, const ScalarFieldd& y) {
            Eigen::VectorXd der(y.size());
            for (int c : interior) {
                der[c] = op.lap(y, c);
            }
            for (int c : boundary) {
                der[c] = 0;
            }
            return der;
        };

        ScalarFieldd T1(domain.size());
        T1[interior] = 1.0;
        T1[boundary] = 0.0;

        auto integrator = integrators::Explicit::Euler().solve(dv_dt, 0.0, time, dt, T1);
        auto stepper = integrator.begin();
        while (stepper) {
            ++stepper;
        }

        return LinfError(stepper.time(), stepper.value());
    }

    double solveABEuler() {
        auto op = shapes.explicitOperators();

        auto dv_dt = [&](double, const ScalarFieldd& y) {
            Eigen::VectorXd der(y.size());
            for (int c : interior) {
                der[c] = op.lap(y, c);
            }
            for (int c : boundary) {
                der[c] = 0;
            }
            return der;
        };

        ScalarFieldd T1(domain.size());
        T1[interior] = 1.0;
        T1[boundary] = 0.0;

        auto integrator = integrators::ExplicitMultistep::AB1().solve(dv_dt, 0.0, time, dt, T1);
        auto stepper = integrator.begin();
        while (stepper) {
            ++stepper;
        }

        return LinfError(stepper.time(), stepper.value());
    }
};

TEST_F(DiffusionExplicitTest, IntegratorsMatch) {
    double err_basic = solveBasic();
    double err_op = solveOperators();
    double err_rk = solveRKEuler();
    double err_ab = solveABEuler();

    double tol = 1e-5;
    EXPECT_NEAR(err_basic, err_op, tol);
    EXPECT_NEAR(err_basic, err_rk, tol);
    EXPECT_NEAR(err_rk, err_ab, tol);
    EXPECT_NEAR(err_ab, err_basic, tol);
    EXPECT_LE(err_basic, 5e-4);
}

}  // namespace mm
