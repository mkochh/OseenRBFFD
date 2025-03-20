#include <medusa/bits/operators/ExplicitVectorOperators.hpp>

#include <medusa/bits/types/Vec.hpp>
#include <medusa/bits/operators/UniformShapeStorage.hpp>
#include <medusa/bits/operators/computeShapes.hpp>
#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/domains/FindClosest.hpp>
#include <medusa/bits/approximations/Monomials_fwd.hpp>
#include <medusa/bits/approximations/RBFBasis.hpp>
#include <medusa/bits/approximations/WeightFunction_fwd.hpp>
#include <medusa/bits/approximations/Gaussian_fwd.hpp>
#include <medusa/bits/approximations/ScaleFunction.hpp>
#include <medusa/bits/approximations/WLS.hpp>
#include <medusa/bits/approximations/JacobiSVDWrapper_fwd.hpp>
#include <medusa/bits/operators/RaggedShapeStorage.hpp>
#include <medusa/bits/operators/ShapeStorage.hpp>
#include <medusa/bits/types/VectorField.hpp>
#include <medusa/bits/operators/ExplicitOperators.hpp>

#include "gtest/gtest.h"

namespace mm {

namespace {
Vec3d a = {1.2, 3.4, -0.7};
Vec3d b = {-0.5, 2.4, 0.8};
Eigen::Matrix3d M;
}

// f(v) = M*v + a + sin(v)
template <int dim>
Vec<double, dim> f(const Vec<double, dim>& v) {
    return M.topLeftCorner<dim, dim>()*v + a.head<dim>() + std::sin(b.head<dim>().dot(v))*v;
}

template <int dim>
Vec<double, dim> df(const Vec<double, dim>& v, int var) {
    double f = b[var]*std::cos(b.head<dim>().dot(v));
    Vec<double, dim> r = M.col(var).head<dim>() + f*v;
    r[var] += std::sin(b.head<dim>().dot(v));
    return r;
}

template <int dim>
Vec<double, dim> ddf(const Vec<double, dim>& v, int varmin, int varmax) {
    double bv = b.head<dim>().dot(v);
    Vec<double, dim> r = -b[varmin]*b[varmax]*std::sin(bv)*v;
    if (varmin == varmax) {
        r[varmin] += 2*b[varmax] * std::cos(bv);
    } else {
        r[varmin] += b[varmax] * std::cos(bv);
        r[varmax] += b[varmin] * std::cos(bv);
    }
    return r;
}

template <int dim>
class ExplicitVecOp : public ::testing::Test {
  public:
    UniformShapeStorage<Vec<double, dim>> storage;
    ExplicitVectorOperators<decltype(storage)> op;
    Range<Vec<double, dim>> pos;
    Range<int> type;
    Range<Vec<double, dim>> field;

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
        typename sh::operator_tuple<sh::all, dim>::type operators{};
        storage.resize(d.supportSizes());
        computeShapes(d, wls, d.all(), operators, &storage);
        op.setShapes(storage);

        // set up test field
        M << 2.4, 1.4, -0.3, 2.8, -1.9, 0.4, -0.3, 2.4, 1.0;
        for (const auto& x : pos) {
            field.push_back(f<dim>(x));
        }
    }
};
typedef ExplicitVecOp<1> ExplicitVecOp1d;
typedef ExplicitVecOp<2> ExplicitVecOp2d;
typedef ExplicitVecOp<3> ExplicitVecOp3d;

TEST_F(ExplicitVecOp1d, Grad) {
    for (int i = 0; i < pos.size(); ++i) {
        Eigen::Matrix<double, 1, 1> grad = op.grad(field, i);
        Vec1d expected = df<1>(pos[i], 0);
        EXPECT_NEAR(expected[0], grad(0, 0), 1e-4);
    }
}

TEST_F(ExplicitVecOp1d, Div) {
    for (int i = 0; i < pos.size(); ++i) {
        if (type[i] < 0) continue;
        double div = op.div(field, i);
        Vec1d expected = df<1>(pos[i], 0);
        EXPECT_NEAR(expected[0], div, 1e-4);
    }
}

TEST_F(ExplicitVecOp1d, Lap) {
    for (int i = 0; i < pos.size(); ++i) {
        if (type[i] < 0) continue;
        Vec1d lap = op.lap(field, i);
        Vec1d expected = ddf<1>(pos[i], 0, 0);
        EXPECT_NEAR(expected[0], lap[0], 1e-4);
    }
}

TEST_F(ExplicitVecOp1d, GradDiv) {
    for (int i = 0; i < pos.size(); ++i) {
        if (type[i] < 0) continue;
        Vec1d gd = op.graddiv(field, i);
        Vec1d expected = ddf<1>(pos[i], 0, 0);
        EXPECT_NEAR(expected[0], gd[0], 1e-4);
    }
}

TEST_F(ExplicitVecOp2d, Grad) {
    for (int i = 0; i < pos.size(); ++i) {
        Eigen::Matrix<double, 2, 2> grad = op.grad(field, i);
        for (int j = 0; j < 2; ++j) {
            Vec2d expected = df<2>(pos[i], j);
            for (int k = 0; k < 2; ++k) {
                EXPECT_NEAR(expected[k], grad(k, j), 1e-2);
            }
        }
    }
}

TEST_F(ExplicitVecOp2d, Div) {
    for (int i = 0; i < pos.size(); ++i) {
        if (type[i] < 0) continue;
        double div = op.div(field, i);
        Vec2d e1 = df<2>(pos[i], 0);
        Vec2d e2 = df<2>(pos[i], 1);
        EXPECT_NEAR(e1[0] + e2[1], div, 2e-3);
    }
}

TEST_F(ExplicitVecOp2d, Lap) {
    for (int i = 0; i < pos.size(); ++i) {
        if (type[i] < 0) continue;
        Vec2d lap = op.lap(field, i);
        Vec2d expected = 0.0;
        for (int k = 0; k < 2; ++k) {
            expected += ddf<2>(pos[i], k, k);
        }
        EXPECT_LT((expected-lap).norm(), 5e-3);
    }
}

TEST_F(ExplicitVecOp2d, GradDiv) {
    for (int i = 0; i < pos.size(); ++i) {
        if (type[i] < 0) continue;
        /// [Explicit syntax]
        Vec2d gd = op.graddiv(field, i);
        /// [Explicit syntax]
        Vec2d expected = 0.0;
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                expected[j] += ddf<2>(pos[i], j, k)[k];
            }
        }
        EXPECT_LT((expected-gd).norm(), 5e-3);
    }
}

TEST_F(ExplicitVecOp3d, Grad) {
    for (int i = 0; i < pos.size(); ++i) {
        Eigen::Matrix<double, 3, 3> grad = op.grad(field, i);
        for (int j = 0; j < 3; ++j) {
            Vec3d expected = df<3>(pos[i], j);
            for (int k = 0; k < 3; ++k) {
                EXPECT_NEAR(expected[k], grad(k, j), 1e-2);
            }
        }
    }
}

TEST_F(ExplicitVecOp3d, Div) {
    for (int i = 0; i < pos.size(); ++i) {
        if (type[i] < 0) continue;
        double div = op.div(field, i);
        Vec3d e1 = df<3>(pos[i], 0);
        Vec3d e2 = df<3>(pos[i], 1);
        Vec3d e3 = df<3>(pos[i], 2);
        EXPECT_NEAR(e1[0] + e2[1] + e3[2], div, 2e-3);
    }
}

TEST_F(ExplicitVecOp3d, Curl) {
    for (int i = 0; i < pos.size(); ++i) {
        if (type[i] <0) continue;
        Vec3d rot = op.curl(field, i);
        Vec3d expected = 0.0;

        Vec3d ux = df<3>(pos[i], 0);
        Vec3d uy = df<3>(pos[i], 1);
        Vec3d uz = df<3>(pos[i], 2);

        expected[0] = uy[2] - uz[1];
        expected[1] = uz[0] - ux[2];
        expected[2] = ux[1] - uy[0];

        for (int k = 0; k < 3; ++k) {
            EXPECT_NEAR(expected[k], rot[k], 1e-2);
        }
    }
}

TEST_F(ExplicitVecOp3d, Lap) {
    for (int i = 0; i < pos.size(); ++i) {
        if (type[i] < 0) continue;
        Vec3d lap = op.lap(field, i);
        Vec3d expected = 0.0;
        for (int k = 0; k < 3; ++k) {
            expected += ddf<3>(pos[i], k, k);
        }
        EXPECT_LT((expected-lap).norm(), 5e-3);
    }
}

TEST_F(ExplicitVecOp3d, GradDiv) {
    for (int i = 0; i < pos.size(); ++i) {
        if (type[i] < 0) continue;
        Vec3d gd = op.graddiv(field, i);
        Vec3d expected = 0.0;
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                expected[j] += ddf<3>(pos[i], j, k)[k];
            }
        }
        EXPECT_LT((expected-gd).norm(), 5e-3);
    }
}

TEST(ExplicitVecOp, GradDifferentDimension) {
    // Tests computing gradient of vector field with different dimension than domain dimension.
    BoxShape<Vec3d> box(0.0, 1.0);
    DomainDiscretization<Vec3d> domain = box.discretizeWithStep(0.25);

    int N = domain.size();
    domain.findSupport(FindClosest(9));

    WLS<Gaussians<Vec3d>, GaussianWeight<Vec3d>, ScaleToFarthest> wls({9, 30.0}, {1.0});
    auto storage = domain.computeShapes(wls);

    VectorField1d vf(N);
    for (int i = 0; i < N; ++i) vf(i, 0) = domain.pos(i, 0);

    auto eop = storage.explicitOperators();
    auto evop = storage.explicitVectorOperators();
    for (int i : domain.interior()) {
        auto egrad = eop.grad(vf.c(0), i);
        auto evgrad = evop.grad(vf, i).transpose();
        EXPECT_EQ(evgrad, egrad);
    }
}


TEST(Operators, CustomExplicitVector) {
    BoxShape<Vec2d> box({0.0, 0.0}, {1.0, 1.0});
    DomainDiscretization<Vec2d> d = box.discretizeWithStep(0.02);
    d.findSupport(FindClosest(9));
    WLS<Monomials<Vec2d>, NoWeight<Vec2d>, ScaleToFarthest> wls(2);

    /// [Custom vector operators]
    struct Dx : Der1<2> { Dx() : Der1(0) {} };
    struct Dy : Der1<2> { Dy() : Der1(1) {} };

    auto storage = d.computeShapes<std::tuple<Dx, Dy, Lap<2>>>(wls);
    ExplicitVectorOperators<decltype(storage)> op(storage);
    VectorField<std::complex<double>, 3> field(d.size());
    std::complex<double> c(2.4, -0.3);
    field.setConstant(c);

    Vec<std::complex<double>, 3> lap = op.apply<Lap<2>>(field, 52);
    Vec<std::complex<double>, 3> dx = op.apply<Dx>(field, 34);
    Vec<std::complex<double>, 3> dy = op.apply<Dy>(field, 32);
    /// [Custom vector operators]
    EXPECT_LT(lap.norm(), 1e-11);
    EXPECT_LT(dx.norm(), 1e-13);
    EXPECT_LT(dy.norm(), 1e-13);
}

/**
 * This test Neumann method on all points in 1D scalar sin field.
 * By by calling the method such that derivation be cos(x), it compares the
 * returned value with sin(x).
 */
TEST(Operators, VectorNeumannSin) {
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
    ExplicitVectorOperators<decltype(storage)> op(storage);

    // Create scalar field of sin(x)
    Range<Vec1d> field(domain.size(), 1);
    for (int i : domain.types() != 0) {
        field[i] = {std::sin(domain.pos(i, 0))};
    }
    for (int i : domain.interior()) {
        double tmp = std::cos(domain.pos(i, 0));
        tmp = op.neumann(field, i, {1}, tmp)[0];
        EXPECT_NEAR(tmp, field[i][0], field[i][0] * 1e-3);
    }
}

TEST(Operators, VectorBoundaryDeathTest) {
    // Create rectangle domain for testing fill it uniformly
    BoxShape<Vec1d> box(0.0, 1.0);
    double dx = 0.2;
    DomainDiscretization<Vec1d> domain = box.discretizeWithStep(dx);

    int supp_size = 3;
    domain.findSupport(FindClosest(supp_size));
    WLS<Monomials<Vec1d>, GaussianWeight<Vec1d>> mls(1, 30*dx);

    auto storage = domain.computeShapes<sh::d1>(mls, domain.interior());
    ExplicitVectorOperators<decltype(storage)> op(storage);

    // Create scalar field
    Range<Vec1d> field(domain.size(), 1);
    for (int i : domain.types() != 0) {
        field[i] = {std::sin(domain.pos(i, 0))};
    }
    for (int i : domain.interior()) {
        EXPECT_DEATH(op.neumann(field, i, {1}, 0.0), "no effect on the flux in direction");
    }
}

TEST(Operators, ExplicitVectorComplexNumbers) {
    BoxShape<Vec2d> box({0.0, 0.0}, {1.0, 1.0});
    DomainDiscretization<Vec2d> d = box.discretizeWithStep(0.02);
    d.findSupport(FindClosest(9));
    WLS<Monomials<Vec2d>, NoWeight<Vec2d>, ScaleToFarthest> wls(2);
    auto storage = d.computeShapes(wls);
    ExplicitVectorOperators<decltype(storage)> op(storage);

    VectorField<std::complex<double>, 2> field(d.size());
    std::complex<double> c1(1.2, -0.4);
    std::complex<double> c2(-1.7, 0.2);
    field = Vec<std::complex<double>, 2>(c1, c2);

    std::complex<double> c(2.4, -0.3);
    field.setConstant(c);

    Vec<std::complex<double>, 2> lap = op.lap(field, 2);
    Vec<std::complex<double>, 2> neu = op.neumann(field, 0, {-1, 3}, std::complex<double>(0.0));
    Eigen::Matrix2cd grad = op.grad(field, 0);
    std::complex<double> div = op.div(field, 0);
    Vec<std::complex<double>, 2> gd = op.graddiv(field, 0);

    EXPECT_LT(lap.norm(), 1e-11);
    EXPECT_LT(gd.norm(), 1e-11);
    EXPECT_LT((neu-Eigen::Vector2cd(c, c)).norm(), 1e-13);
    EXPECT_LT(grad.norm(), 1e-13);
    EXPECT_LT(std::abs(div), 1e-13);
}

TEST(Operators, DISABLED_ExplicitVectorOperatorsUsageExample) {
    /// [Explicit vector operators usage example]
    BoxShape<Vec2d> box({0.0, 0.0}, {1.0, 1.0});
    DomainDiscretization<Vec2d> d = box.discretizeWithStep(0.02);
    d.findSupport(FindClosest(9));
    WLS<Monomials<Vec2d>, NoWeight<Vec2d>, ScaleToFarthest> wls(2);
    auto storage = d.computeShapes(wls);
    ExplicitVectorOperators<decltype(storage)> op(storage);
    std::cout << op << std::endl;
    VectorField2d v(d.size());  // constant scalar field, value 1.34 at each point
    v = Vec2d(12.4, -3.4);
    Vec2d val = op.lap(v, 2);  // laplace at point d.pos(2);
    Eigen::Matrix2d grad = op.grad(v, 4);  // gradient at point d.pos(4);
    double div = op.div(v, 3);  // divergence at points d.pos(3);

    VectorField<std::complex<double>, 3> complex_field(d.size());  // 3d complex field
    complex_field.setConstant(std::complex<double>(2.5, -7.1));
    Eigen::Matrix<std::complex<double>, 3, 2> cgrad = op.grad(complex_field, 1);
    /// [Explicit vector operators usage example]
    (void) div; (void) val; (void) cgrad; (void) grad;
}

}  // namespace mm
