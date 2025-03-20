#include <medusa/bits/approximations/Monomials_fwd.hpp>
#include <medusa/bits/approximations/Gaussian_fwd.hpp>
#include <medusa/bits/operators/RaggedShapeStorage.hpp>
#include <medusa/bits/approximations/RBFFD_fwd.hpp>
#include <medusa/bits/io/HDF.hpp>
#include <medusa/bits/io/HDF_Eigen.hpp>
#include <medusa/bits/domains/DomainDiscretization_fwd.hpp>
#include <medusa/bits/domains/BoxShape_fwd.hpp>
#include <medusa/bits/domains/FindClosest.hpp>
#include <medusa/bits/operators/UniformShapeStorage_fwd.hpp>
#include <medusa/bits/operators/computeShapes.hpp>
#include <medusa/bits/operators/ImplicitVectorOperators.hpp>
#include <medusa/bits/operators/ExplicitVectorOperators.hpp>
#include <medusa/bits/types/ScalarField.hpp>
#include <medusa/bits/types/VectorField.hpp>
#include <medusa/bits/approximations/WeightFunction_fwd.hpp>
#include <medusa/bits/approximations/ScaleFunction.hpp>
#include <medusa/bits/approximations/WLS.hpp>
#include <medusa/bits/approximations/JacobiSVDWrapper_fwd.hpp>
#include "Eigen/SparseCore"
#include "Eigen/IterativeLinearSolvers"

#include "gtest/gtest.h"

namespace mm {

class CantileverBeamTest : public ::testing::Test {
  public:
    // Physical parameters
    static constexpr double E = 72.1e9;
    static constexpr double v = 0.33;
    static constexpr double P = 1000;
    static constexpr double D = 5;
    static constexpr double L = 30;

    // Numerical parameters
    static constexpr double solver_tolerance = 1e-15;
    static constexpr int max_iterations = 50;
    static constexpr double drop_tolerance = 1e-5;
    static constexpr int fill_factor = 20;

    typedef Eigen::SparseMatrix<double, Eigen::RowMajor> matrix_t;
    Eigen::BiCGSTAB<matrix_t, Eigen::IncompleteLUT<double>> solver;

    template <class BasisType>
    double solve(const BasisType& basis, int num_nodes_x, double sigmaW = 1.0) {
        // Derived parameters
        const double I = D * D * D / 12;
        double lam = E * v / (1 - 2 * v) / (1 + v);
        const double mu = E / 2 / (1 + v);
        lam = 2 * mu * lam / (2 * mu + lam);  // plane stress

        // Closed form solution -- refer to cantilever_beam.nb for reference.
        std::function<Vec2d(Vec2d)> analytical = [=](const Vec2d& p) {
            double x = p[0], y = p[1];
            double ux = (P*y*(3*D*D*(1+v)-4*(3*L*L-3*x*x+(2+v)*y*y))) / (2.*D*D*D*E);
            double uy = -(P*(3*D*D*(1+v)*(L-x) + 4*(L-x)*(L-x)*(2*L+x)+12*v*x*y*y)) / (2.*D*D*D*E);
            return Vec2d(ux, uy);
        };

        int BOTTOM = -3;
        int TOP = -4;
        int LEFT = -1;
        int RIGHT = -2;

        // Domain definition
        Vec2d low(0, -D / 2), high(L, D / 2);
        BoxShape<Vec2d> box(low, high);
        double step = L / num_nodes_x;
        DomainDiscretization<Vec2d> domain = box.discretizeWithStep(step);
        int N = domain.size();

        int support_size = 9;
        domain.findSupport(FindClosest(support_size));

        WLS<BasisType, GaussianWeight<Vec2d>, ScaleToClosest> approx(basis, sigmaW);
        auto storage = domain.computeShapes(approx);

        Eigen::SparseMatrix<double, Eigen::RowMajor> M(2*N, 2*N);
        M.reserve(storage.supportSizesVec());

        Eigen::VectorXd rhs(2*N); rhs.setZero(); rhs.setZero();
        auto op = storage.implicitVectorOperators(M, rhs);

        // Set equation on interior
        for (int i : domain.interior()) {
            (lam+mu)*op.graddiv(i) + mu*op.lap(i) = 0.0;
        }
        for (int i : domain.types() == RIGHT) {
            double y = domain.pos(i, 1);
            op.value(i) = {(P*y*(3*D*D*(1+v) - 4*(2+v)*y*y)) / (24.*E*I), -(L*v*P*y*y) / (2.*E*I)};
        }
        for (int i : domain.types() == LEFT) {
            double y = domain.pos(i, 1);
            op.traction(i, lam, mu, {-1, 0}) = {0, -P*(D*D - 4*y*y) / (8.*I)};
        }
        for (int i : domain.types() == TOP) {
            op.traction(i, lam, mu, {0, 1}) = 0.0;
        }
        for (int i : domain.types() == BOTTOM) {
            op.traction(i, lam, mu, {0, -1}) = 0.0;
        }
        M.makeCompressed();

        solver.preconditioner().setDroptol(drop_tolerance);
        solver.preconditioner().setFillfactor(fill_factor);
        solver.setMaxIterations(max_iterations);
        solver.setTolerance(solver_tolerance);
        solver.compute(M);
        Eigen::VectorXd sol = solver.solve(rhs);

        Range<double> error(N);
        double maxuv = 0;
        for (int i = 0; i < N; ++i) {
            Vec2d uv = analytical(domain.pos(i));
            maxuv = std::max(maxuv, std::max(std::abs(uv[0]), std::abs(uv[1])));
            error[i] = std::max(std::abs(uv[0] - sol[i]), std::abs(uv[1] - sol[i + N]));
        }

        double L_inf_error = *std::max_element(error.begin(), error.end()) / maxuv;
        return L_inf_error;
    }
};
constexpr double CantileverBeamTest::solver_tolerance;
constexpr double CantileverBeamTest::drop_tolerance;
constexpr double CantileverBeamTest::L;

TEST_F(CantileverBeamTest, Monomials) {
    Monomials<Vec2d> mon9 = Monomials<Vec2d>::tensorBasis(2);

    int nx = 100;
    double err = solve(mon9, nx);
    EXPECT_LE(err, 1.9e-2);
    EXPECT_LE(solver.iterations(), 8);
    EXPECT_LE(solver.error(), solver_tolerance);

    nx = 200;
    err = solve(mon9, nx);
    EXPECT_LE(err, 5e-3);
    EXPECT_LE(solver.iterations(), 12);
    EXPECT_LE(solver.error(), solver_tolerance);
}

TEST_F(CantileverBeamTest, MonomialsNoWeight) {
    Monomials<Vec2d> mon9 = Monomials<Vec2d>::tensorBasis(2);
    int nx = 100;  /* unstable at first, only test later */
    double err = solve(mon9, nx, 10000.0);
    EXPECT_LE(err, 2.5e-2);
    EXPECT_LE(solver.iterations(), 8);
    EXPECT_LE(solver.error(), solver_tolerance);

    nx = 200;
    err = solve(mon9, nx, 10000.0);
    EXPECT_LE(err, 5e-3);
    EXPECT_LE(solver.iterations(), 20);
    EXPECT_LE(solver.error(), solver_tolerance);
}

TEST_F(CantileverBeamTest, Gaussians) {
    int nx = 100;
    double sigmaB = 400.0;
    Gaussians<Vec2d> gau9(9, sigmaB);
    double err = solve(gau9, nx);
    EXPECT_LE(err, 2.5e-2);
    EXPECT_LE(solver.iterations(), 8);
    EXPECT_LE(solver.error(), solver_tolerance);

    nx = 200;
    err = solve(gau9, nx);
    EXPECT_LE(err, 8e-3);
    EXPECT_LE(solver.iterations(), 12);
    EXPECT_LE(solver.error(), solver_tolerance);
}

}  // namespace mm
