#include <medusa/bits/operators/ExplicitOperators.hpp>
#include <medusa/bits/types/Range_fwd.hpp>
#include <medusa/bits/types/Vec_fwd.hpp>
#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/approximations/WLS.hpp>
#include <medusa/bits/approximations/Monomials.hpp>
#include <medusa/bits/approximations/ScaleFunction.hpp>
#include <medusa/bits/approximations/RBFBasis.hpp>
#include <medusa/bits/approximations/WeightFunction_fwd.hpp>
#include <medusa/bits/approximations/JacobiSVDWrapper_fwd.hpp>
#include <medusa/bits/operators/UniformShapeStorage.hpp>
#include <medusa/bits/operators/RaggedShapeStorage.hpp>
#include <medusa/bits/operators/ShapeStorage.hpp>
#include <medusa/bits/operators/computeShapes.hpp>
#include <medusa/bits/domains/FindClosest.hpp>
#include <medusa/bits/types/ScalarField_fwd.hpp>
#include <medusa/bits/approximations/Gaussian.hpp>

#include <cmath>

#include "gtest/gtest.h"

namespace mm {

namespace {
Vec3d a = {1.2, 3.4, -0.7};
Vec3d b = {-0.8, 2.3, 1.4};
double c = -0.4;
}

// f(x, y, z) = sin(b.x) + a.x + c
template <int dim>
double f(const Vec<double, dim>& v) {
    return std::sin(b.head<dim>().dot(v)) + a.head<dim>().dot(v) + c;
}

template <int dim>
double df(const Vec<double, dim>& v, int var) {
    return b[var] * std::cos(b.head<dim>().dot(v)) + a[var];
}

template <int dim>
double ddf(const Vec<double, dim>& v, int var1, int var2) {
    return -b[var1]*b[var2]*std::sin(b.head<dim>().dot(v));
}

template <int dim>
class ExplicitOp : public ::testing::Test {
  public:
    UniformShapeStorage<Vec<double, dim>> storage;
    ExplicitOperators<decltype(storage)> op;
    Range<Vec<double, dim>> pos;
    Range<int> type;
    Range<double> field;

  protected:
    void SetUp() override {
        // set up operators
        BoxShape<Vec<double, dim>> box(2.23, 2.34);
        DomainDiscretization<Vec<double, dim>> d = box.discretizeWithStep(0.02);
        pos = d.positions();
        type = d.types();
        int num_cl = 1; for (int i = 0; i < dim; ++i) num_cl *= 3;  // num_cl = 3^dim
        d.findSupport(FindClosest(num_cl));
        WLS<Monomials<Vec<double, dim>>, NoWeight<Vec<double, dim>>, ScaleToFarthest>
        wls(Monomials<Vec<double, dim>>::tensorBasis(2));
        typename sh::operator_tuple<sh::all, dim>::type operators;
        storage.resize(d.supportSizes());
        computeShapes(d, wls, d.all(), operators, &storage);
        op.setShapes(storage);

        // set up test field
        for (const auto& x : pos) {
            field.push_back(f<dim>(x));
        }
    }
};
typedef ExplicitOp<1> ExplicitOp1d;
typedef ExplicitOp<2> ExplicitOp2d;
typedef ExplicitOp<3> ExplicitOp3d;

TEST_F(ExplicitOp1d, Grad) {
    for (int i = 0; i < pos.size(); ++i) {
        Vec1d grad = op.grad(field, i);
        double expected = df<1>(pos[i], 0);
        EXPECT_NEAR(expected, grad[0], 1e-4);
    }
}

TEST_F(ExplicitOp1d, Lap) {
    for (int i = 0; i < pos.size(); ++i) {
        if (type[i] < 0) continue;
        double lap = op.lap(field, i);
        double expected = ddf<1>(pos[i], 0, 0);
        EXPECT_NEAR(expected, lap, 1e-4);
    }
}

TEST_F(ExplicitOp1d, D1) {
    for (int i = 0; i < pos.size(); ++i) {
        double der = op.d1(field, 0, i);
        double expected = df<1>(pos[i], 0);
        EXPECT_NEAR(expected, der, 1e-4);
    }
}

TEST_F(ExplicitOp1d, D2) {
    for (int i = 0; i < pos.size(); ++i) {
        if (type[i] < 0) continue;
        double der = op.d2(field, 0, 0, i);
        double expected = ddf<1>(pos[i], 0, 0);
        EXPECT_NEAR(expected, der, 1e-4);
    }
}

TEST_F(ExplicitOp2d, Grad) {
    for (int i = 0; i < pos.size(); ++i) {
        Vec2d grad = op.grad(field, i);
        for (int j = 0; j < 2; ++j) {
            double expected = df<2>(pos[i], j);
            EXPECT_NEAR(expected, grad[j], 5e-3);
        }
    }
}

TEST_F(ExplicitOp2d, Lap) {
    for (int i = 0; i < pos.size(); ++i) {
        if (type[i] < 0) continue;
        double lap = op.lap(field, i);
        double lap_exp = 0;
        for (int j = 0; j < 2; ++j) {
            lap_exp += ddf<2>(pos[i], j, j);
        }
        EXPECT_NEAR(lap_exp, lap, 5e-3);
    }
}

TEST_F(ExplicitOp2d, D1) {
    for (int d = 0; d < 2; ++d) {
        for (int i = 0; i < pos.size(); ++i) {
            double der = op.d1(field, d, i);
            double expected = df<2>(pos[i], d);
            EXPECT_NEAR(expected, der, 5e-3);
        }
    }
}

TEST_F(ExplicitOp2d, D2) {
    for (int dmax = 0; dmax < 2; ++dmax) {
        for (int dmin = 0; dmin <= dmax; ++dmin) {
            for (int i = 0; i < pos.size(); ++i) {
                if (type[i] < 0) continue;
                double der = op.d2(field, dmin, dmax, i);
                double expected = ddf<2>(pos[i], dmin, dmax);
                EXPECT_NEAR(expected, der, 5e-3);
            }
        }
    }
}

TEST_F(ExplicitOp3d, Grad) {
    for (int i = 0; i < pos.size(); ++i) {
        Vec3d grad = op.grad(field, i);
        for (int j = 0; j < 3; ++j) {
            double expected = df<3>(pos[i], j);
            EXPECT_NEAR(expected, grad[j], 5e-3);
        }
    }
}

TEST_F(ExplicitOp3d, Lap) {
    for (int i = 0; i < pos.size(); ++i) {
        if (type[i] < 0) continue;
        double lap = op.lap(field, i);
        double lap_exp = 0;
        for (int j = 0; j < 3; ++j) {
            lap_exp += ddf<3>(pos[i], j, j);
        }
        EXPECT_NEAR(lap_exp, lap, 1e-3);
    }
}

TEST_F(ExplicitOp3d, D1) {
    for (int d = 0; d < 3; ++d) {
        for (int i = 0; i < pos.size(); ++i) {
            double der = op.d1(field, d, i);
            double expected = df<3>(pos[i], d);
            EXPECT_NEAR(expected, der, 5e-3);
        }
    }
}

TEST_F(ExplicitOp3d, D2) {
    for (int dmax = 0; dmax < 3; ++dmax) {
        for (int dmin = 0; dmin <= dmax; ++dmin) {
            for (int i = 0; i < pos.size(); ++i) {
                if (type[i] < 0) continue;
                double der = op.d2(field, dmin, dmax, i);
                double expected = ddf<3>(pos[i], dmin, dmax);
                EXPECT_NEAR(expected, der, 1e-3);
            }
        }
    }
}

/**
 * This test Neumann method on all points in 1D scalar sin field.
 * By by calling the method such that derivation be cos(x), it compares the
 * returned value with sin(x).
 */
TEST(Operators, ScalarNeumannSin) {
    // Create rectangle domain for testing fill it uniformly
    BoxShape<Vec1d> box(0.0, 1.0);
    double dx = 0.01;
    DomainDiscretization<Vec1d> domain = box.discretizeWithStep(dx);
    // You need asymmetry in selecting the support so that the derivation is
    // is dependant on the central point
    int supp_size = 4;
    domain.findSupport(FindClosest(supp_size));
    WLS<Monomials<Vec1d>, GaussianWeight<Vec1d>> mls(1, 30*dx);
    auto storage = domain.computeShapes<sh::d1>(mls, domain.interior());
    ExplicitOperators<decltype(storage)> op(storage);

    // Create scalar field of sin(x)
    Range<double> field(domain.size(), 1);
    for (int i : domain.types() != 0) {
        field[i] = std::sin(domain.pos(i, 0));
    }
    for (int i : domain.interior()) {
        double tmp = std::cos(domain.pos(i, 0));
        tmp = op.neumann(field, i, {1}, tmp);
        EXPECT_NEAR(tmp, field[i], field[i] * 1e-3);
    }
}

TEST(Operators, ScalarBoundaryDeathTest) {
    // Create rectangle domain for testing fill it uniformly
    BoxShape<Vec1d> box(0.0, 1.0);
    double dx = 0.2;
    DomainDiscretization<Vec1d> domain = box.discretizeWithStep(dx);

    int supp_size = 3;
    domain.findSupport(FindClosest(supp_size));
    WLS<Monomials<Vec1d>, GaussianWeight<Vec1d>> mls(1, 30*dx);

    auto storage = domain.computeShapes<sh::d1>(mls, domain.interior());
    ExplicitOperators<decltype(storage)> op(storage);

    // Create scalar field
    Range<double> field(domain.size(), 1);
    for (int i : domain.types() != 0) {
        field[i] = std::sin(domain.pos(i, 0));
    }
    for (int i : domain.interior()) {
        EXPECT_DEATH(op.neumann(field, i, {1}, 0.0), "no effect on the flux in direction");
    }
}

/**
 * @brief Closed form solution of diffusion equation
 * @param pos spatial coordinate
 * @param t time
 * @param a size of domain
 * @param D diffusion constant
 * @param N no. of expansion
 * @return value of temperature
 */
double diffusion_analytical(const Vec2d& pos, double t, double a, double D, size_t N) {
    double T = 0;
    double f = PI / a;
    for (size_t n = 1; n < N; n = n + 2) {
        for (size_t m = 1; m < N; m = m + 2) {
            T += 4.0 / f / f / (n * m) * std::sin(n * f * pos[0]) * std::sin(m * f * pos[1]) *
                 std::exp(-D * t * ((n * n + m * m) * f * f));
        }
    }
    return T;
}

// A test for Neumann boundary conditions
TEST(Operators, QuarterDiffusionTest) {
    // Create rectangle domain for testing fill it uniformly
    BoxShape<Vec2d> box(0.0, 1.0);
    double dx = 0.05;
    DomainDiscretization<Vec2d> domain = box.discretizeWithStep(dx);

    int N = domain.size();
    auto internal = domain.interior();
    auto up_right = domain.types().filter([](int i) { return i == -2 || i == -4; });
    auto down_left = domain.types().filter([](int i) { return i == -1 || i == -3; });
    auto testing = domain.positions().filter([](const Vec2d& p) {
        return p[0] != 0.0 && p[1] != 1.0 && p[0] < 0.5 && p[1] < 0.5; });

    int supp_size = 5;
    // Find normal support for internal nodes
    FindClosest f(supp_size); f.forNodes(internal);
    domain.findSupport(f);

    f.forNodes(down_left).searchAmong(internal+up_right).forceSelf();
    domain.findSupport(f);

    // Make MLS with Gaussian basis and weight
    WLS<Gaussians<Vec2d>, GaussianWeight<Vec2d>> wls({5, 10*dx}, 30*dx);

    auto storage = domain.computeShapes<sh::lap|sh::d1>(wls, internal+down_left);
    ExplicitOperators<decltype(storage)> op(storage);

    // Create scalar field
    ScalarFieldd u1(N), u2(N);
    u1 = 1.0;
    u2 = 1.0;
    u1[up_right] = 0.0;

    double t = 0.1, dt = 1e-5;
    for (int nt = 0; nt*dt < t; nt++) {
        // Test the correctness of the solution
        if (nt % static_cast<int>(t / dt / 100) == 0 && nt != 0) {
            double tmp_val;
            for (int i : testing) {
                tmp_val = diffusion_analytical(
                        domain.pos(i) + Vec2d(1, 1), nt * dt, 2, 1, 50);
                EXPECT_NEAR(tmp_val, u1[i], 6e-3);
            }
        }
        // Calculate internal nodes with Laplace operator
        for (int i : internal) {
            u2[i] += dt * op.lap(u1, i);
        }
        u1[internal] = u2[internal];
        // Calculate Neumann boundary condition
        for (int i : down_left) {
            u1[i] = op.neumann(u1, i, domain.normal(i), 0.0);
        }
    }
}

TEST(Operators, ExplicitComplexNumbers) {
    BoxShape<Vec2d> box({0.0, 0.0}, {1.0, 1.0});
    DomainDiscretization<Vec2d> d = box.discretizeWithStep(0.02);
    d.findSupport(FindClosest(9));
    WLS<Monomials<Vec2d>, NoWeight<Vec2d>, ScaleToFarthest> wls(2);
    auto storage = d.computeShapes(wls);
    ExplicitOperators<decltype(storage)> op(storage);
    ScalarField<std::complex<double>> field(d.size());
    std::complex<double> c(2.4, -0.3);
    field.setConstant(c);

    std::complex<double> lap = op.lap(field, 2);
    std::complex<double> neu = op.neumann(field, 0, {-1}, 0.0 /* this 0.0 is cast to complex */);
    EXPECT_LT(std::abs(lap), 1e-11);
    EXPECT_LT(std::abs(neu - c), 1e-13);
}

TEST(Operators, CustomExplicit) {
    BoxShape<Vec2d> box({0.0, 0.0}, {1.0, 1.0});
    DomainDiscretization<Vec2d> d = box.discretizeWithStep(0.02);
    d.findSupport(FindClosest(9));
    WLS<Monomials<Vec2d>, NoWeight<Vec2d>, ScaleToFarthest> wls(2);

    /// [Custom explicit operators]
    struct Dx : Der1<2> { Dx() : Der1(0) {} };
    struct Dy : Der1<2> { Dy() : Der1(1) {} };

    auto storage = d.computeShapes<std::tuple<Dx, Dy, Lap<2>>>(wls);
    ExplicitOperators<decltype(storage)> op(storage);
    ScalarField<std::complex<double>> field(d.size());
    std::complex<double> c(2.4, -0.3);
    field.setConstant(c);

    std::complex<double> lap = op.apply<Lap<2>>(field, 52);
    std::complex<double> dx = op.apply<Dx>(field, 34);
    std::complex<double> dy = op.apply<Dy>(field, 32);
    /// [Custom explicit operators]
    EXPECT_LT(std::abs(lap), 1e-11);
    EXPECT_LT(std::abs(dx), 1e-13);
    EXPECT_LT(std::abs(dy), 1e-13);
}

TEST(Operators, DISABLED_ExplicitUsageExample) {
    /// [Explicit operators usage example]
    BoxShape<Vec2d> box({0.0, 0.0}, {1.0, 1.0});
    DomainDiscretization<Vec2d> d = box.discretizeWithStep(0.02);
    d.findSupport(FindClosest(9));
    WLS<Monomials<Vec2d>, NoWeight<Vec2d>, ScaleToFarthest> wls(2);
    auto storage = d.computeShapes(wls);
    ExplicitOperators<decltype(storage)> op(storage);
    std::cout << op << std::endl;
    Range<double> field(d.size(), 1.34);  // constant scalar field, value 1.34 at each point
    double val = op.lap(field, 2);  // laplace at point d.pos(2);
    Vec2d grad = op.grad(field, 4);  // gradient at point d.pos(4);

    ScalarField<std::complex<double>> complex_field(d.size());
    std::complex<double> c(2.4, -0.3);
    complex_field.setConstant(c);  // constant complex field, value c at each point
    // what should the value at node 0 be, so that the normal derivative is 0 + 0i?
    std::complex<double> neu = op.neumann(field, 0, {-1}, 0.0);  // result equal to c
    /// [Explicit operators usage example]
    (void) grad;
    (void) val;
    (void) neu;
}

}  // namespace mm
