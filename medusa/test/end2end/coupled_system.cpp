#include <medusa/Medusa_fwd.hpp>
#include <complex>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

#include "gtest/gtest.h"

namespace mm {

TEST(End2end, CoupledSystem) {
    BoxShape<Vec2d> box(0.0, 1.0);
    DomainDiscretization<Vec2d> domain = box.discretizeWithStep(0.1);

    int N = domain.size();
    domain.findSupport(FindClosest(9));

    WLS<Gaussians<Vec2d>> wls({9, 30.0}, {});
    auto storage = domain.computeShapes<0>(wls);

    // Implicit scalar
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(2*N, 2*N);
    M.reserve(Range<int>(2*N, 2));
    Eigen::VectorXd rhs(2*N); rhs.setZero();

    auto op1 = storage.implicitOperators(M, rhs);
    auto op2 = storage.implicitOperators(M, rhs);
    op2.setRowOffset(N);
    op2.setColOffset(N);
    for (int i : domain.all()) {
        op1.value(i) + op2.value(i, i-N) = 0;
        op1.value(i, i+N) + -op2.value(i) = 2;
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.compute(M);
    Eigen::VectorXd u = solver.solve(rhs);

    EXPECT_NEAR((u.head(N) - Eigen::VectorXd::Ones(N)).norm(), 0, 0);
    EXPECT_NEAR((u.tail(N) + Eigen::VectorXd::Ones(N)).norm(), 0, 0);
}

TEST(End2end, CoupledVectorSystem) {
    constexpr int dim = 1;
    typedef Vec<double, dim> vec;
    BoxShape<vec> box(0.0, 1.0);
    DomainDiscretization<vec> domain = box.discretizeWithStep(0.1);

    int N = domain.size();
    domain.findSupport(FindClosest(9));

    WLS<Monomials<vec>, GaussianWeight<vec>, ScaleToFarthest> wls(2, 1.0);  // irrelevant
    auto storage = domain.computeShapes<0>(wls);

    // Implicit scalar
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(2*dim*N, 2*dim*N);
    M.reserve(Range<int>(2*dim*N, 2));
    Eigen::VectorXd rhs(2*dim*N); rhs.setZero();

    auto op1 = storage.implicitVectorOperators(M, rhs);
    auto op2 = storage.implicitVectorOperators(M, rhs);
    op2.setRowOffset(dim*N);
    op2.setColOffset(dim*N);
    for (int i : domain.all()) {
        op1.value(i) + op2.value(i, i-dim*N) = 0;
        op1.value(i, i+dim*N) + -op2.value(i) = 2;
    }

    Eigen::BiCGSTAB<decltype(M), Eigen::IncompleteLUT<double>> solver;
    solver.compute(M);
    Eigen::VectorXd u = solver.solve(rhs);

    EXPECT_NEAR((u.head(dim*N) - Eigen::VectorXd::Ones(dim*N)).norm(), 0, 0);
    EXPECT_NEAR((u.tail(dim*N) + Eigen::VectorXd::Ones(dim*N)).norm(), 0, 0);
}

}   // namespace mm
