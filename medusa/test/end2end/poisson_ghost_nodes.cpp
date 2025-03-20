#include <medusa/Medusa_fwd.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/LU>

#include "gtest/gtest.h"

using namespace mm;  // NOLINT(*)
using namespace std;  // NOLINT(*)

TEST(End2end, GhostNodes) {
    /// [Ghost nodes]
    // 3D Poisson problem with mixed BC
    double step = 0.1;
    BallShape<Vec3d> sh(0, 1);
    DomainDiscretization<Vec3d> domain = sh.discretizeBoundaryWithStep(step);

    auto dir = domain.positions().filter([](const Vec3d& p) { return p[0] > 0; });
    auto neu = domain.positions().filter([](const Vec3d& p) { return p[0] <= 0; });

    GeneralFill<Vec3d> fill; fill.seed(20);
    domain.fill(fill, step);

    // Add ghost nodes and remember the mapping from boundary to ghost nodes.
    Range<int> gh = domain.addGhostNodes(step, 0, neu);
    int N = domain.size();

    // By default, support is found for non-zero nodes, searching among all nodes
    int n = 25;
    domain.findSupport(FindClosest(n));
    RBFFD<Polyharmonic<double, 3>, Vec3d, ScaleToClosest> approx({}, 2);
    // Approximations are computed for non-zero nodes.
    auto storage = domain.computeShapes(approx);

    // Matrix size is equal to the count of all nodes.
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(N, N);
    // Warning: support size vector contains zeros in ghost nodes, but we want matrix to allocate
    // storage for ghost nodes as well.
    M.reserve(Range<int>(N, n));

    Eigen::VectorXd rhs(N); rhs.setZero();
    auto op = storage.implicitOperators(M, rhs);

    for (int i : domain.interior()) {
        op.lap(i) = 1;
    }
    for (int i : dir) {
        op.value(i) = 0;
    }
    for (int i : neu) {
        op.neumann(i, domain.normal(i)) = 0;
        op.lap(i, gh[i]) = 1;
    }

    Eigen::SparseLU<decltype(M)> solver;
    solver.compute(M);
    ScalarFieldd u = solver.solve(rhs);
    /// [Ghost nodes]

    EXPECT_NEAR(-0.294074915331889, u.mean(), 1e-11);
    EXPECT_NEAR(0, u.maxCoeff(), 1e-11);
    EXPECT_NEAR(-0.612629726851429, u.minCoeff(), 1e-11);
}
