#include <medusa/bits/types/Vec.hpp>
#include <medusa/bits/utils/numutils.hpp>
#include <medusa/bits/domains/FindClosest.hpp>
#include <medusa/bits/approximations/Polyharmonic.hpp>
#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/approximations/RBFFD.hpp>
#include <Eigen/Core>

#include "gtest/gtest.h"
#include "medusa/bits/interpolants/PUApproximant_fwd.hpp"

namespace mm {

TEST(Interpolants, PUApproximant1D) {
    /// [Partition-of-unity approximant usage example]
    // Domain.
    BoxShape<Vec1d> b(0, 1);
    DomainDiscretization<Vec1d> domain = b.discretizeBoundaryWithStep(0.01);
    int support_size = 2;
    domain.findSupport(FindClosest(support_size));

    // Values.
    Eigen::VectorXd values(2);
    for (int i = 0; i < values.size(); i++) {
        values[i] = static_cast<double>(i);
    }

    // Approximation engine.
    RBFFD<Polyharmonic<double>, Vec1d, ScaleToClosest> approx({3}, Monomials<Vec1d>(2));

    // Test points.
    Range<Vec1d> test_points(1);
    test_points[0] = Vec1d(0.5);

    // Partition of unity approximation.
    double effective_radius = 1.0;  // Effective radius.
    // Interpolated value.
    auto value = PUApproximant<Vec1d>::evaluate(domain, values, test_points, effective_radius,
                                                approx);  // value = 0.5

    /// [Partition-of-unity approximant usage example]
    EXPECT_NEAR(0.5, value[0], 1e-15);
}

TEST(Interpolants, PUApproximant2D) {
    // Domain.
    BoxShape<Vec2d> b(0, 1);
    DomainDiscretization<Vec2d> domain = b.discretizeWithStep(0.05);
    int support_size = 15;
    domain.findSupport(FindClosest(support_size));

    // Values.
    Eigen::VectorXd values(domain.size());
    values.setOnes();

    // Approximation engine.
    RBFFD<Polyharmonic<double>, Vec2d, ScaleToClosest> approx({3}, Monomials<Vec2d>(2));

    // Test points.
    BoxShape<Vec2d> b_test(0, 1);
    DomainDiscretization<Vec2d> domain_test = b_test.discretizeWithStep(0.1);
    auto test_points = domain_test.positions();

    // Partition of unity approximation.
    double effective_radius = 1.0;  // Effective radius.
    // Interpolated value.
    auto interpolated_values =
        PUApproximant<Vec2d>::evaluate(domain, values, test_points, effective_radius, approx);

    for (auto value : interpolated_values) {
        EXPECT_NEAR(1, value, 1e-15);
    }
}

TEST(Interpolants, PUApproximant2DNotConstant) {
    // Domain.
    BoxShape<Vec2d> b(0, 1);
    DomainDiscretization<Vec2d> domain = b.discretizeBoundaryWithStep(1);
    int support_size = 4;
    domain.findSupport(FindClosest(support_size));

    // Values.
    int N = domain.size();
    Eigen::VectorXd values(N);
    for (int i = 0; i < N; i++) {
        values[i] = i;
    }

    // Approximation engine.
    RBFFD<Polyharmonic<double>, Vec2d, ScaleToClosest> approx({3}, Monomials<Vec2d>(2));

    // Test points.
    Range<Vec2d> test_points;
    test_points.push_back(Vec2d(0.5, 0.5));

    // Partition of unity approximation.
    double effective_radius = 1.0;  // Effective radius.
    // Interpolated value.
    auto num_values =
        PUApproximant<Vec2d>::evaluate(domain, values, test_points, effective_radius, approx);

    EXPECT_NEAR(1, num_values[0], 1e-15);
}

TEST(Interpolants, PUApproximant3D) {
    // Domain.
    BoxShape<Vec3d> b(0, 1);
    DomainDiscretization<Vec3d> domain = b.discretizeWithStep(0.1);
    int support_size = 30;
    domain.findSupport(FindClosest(support_size));

    // Values.
    Eigen::VectorXd values(domain.size());
    values.setOnes();

    // Approximation engine.
    RBFFD<Polyharmonic<double>, Vec3d, ScaleToClosest> approx({3}, Monomials<Vec3d>(2));

    // Test points.
    BoxShape<Vec3d> b_test(0, 1);
    DomainDiscretization<Vec3d> domain_test = b_test.discretizeWithStep(0.2);
    auto test_points = domain_test.positions();

    // Partition of unity approximation.
    double effective_radius = 1.0;  // Effective radius.
    // Interpolated value.
    auto interpolated_values =
        PUApproximant<Vec3d>::evaluate(domain, values, test_points, effective_radius, approx);

    for (auto value : interpolated_values) {
        EXPECT_NEAR(1, value, 1e-15);
    }
}

TEST(Interpolants, PUApproximantKDTree1D) {
    /// [Partition-of-unity approx tree usage example]
    // Domain.
    BoxShape<Vec1d> b(0, 1);
    DomainDiscretization<Vec1d> domain = b.discretizeBoundaryWithStep(0.01);
    int support_size = 2;
    domain.findSupport(FindClosest(support_size));

    // Values.
    Eigen::VectorXd values(2);
    for (int i = 0; i < values.size(); i++) {
        values[i] = static_cast<double>(i);
    }

    // Approximation engine.
    RBFFD<Polyharmonic<double>, Vec1d, ScaleToClosest> approx({3}, Monomials<Vec1d>(2));

    // Test points.
    Range<Vec1d> test_points(1);
    test_points[0] = Vec1d(0.5);
    KDTree<Vec1d> tree(test_points);

    // Partition of unity approximation.
    double effective_radius = 1.0;  // Effective radius.
    // Interpolated value.
    auto value = PUApproximant<Vec1d>::evaluate(domain, values, tree, effective_radius,
                                                approx);  // value = 0.5

    /// [Partition-of-unity approx tree usage example]
    EXPECT_NEAR(0.5, value[0], 1e-15);
}

TEST(Interpolants, PUApproximantKDTree2D) {
    // Domain.
    BoxShape<Vec2d> b(0, 1);
    DomainDiscretization<Vec2d> domain = b.discretizeWithStep(0.05);
    int support_size = 15;
    domain.findSupport(FindClosest(support_size));

    // Values.
    Eigen::VectorXd values(domain.size());
    values.setOnes();

    // Approximation engine.
    RBFFD<Polyharmonic<double>, Vec2d, ScaleToClosest> approx({3}, Monomials<Vec2d>(2));

    // Test points.
    BoxShape<Vec2d> b_test(0, 1);
    DomainDiscretization<Vec2d> domain_test = b_test.discretizeWithStep(0.1);
    auto test_points = domain_test.positions();
    KDTree<Vec2d> tree(test_points);

    // Partition of unity approximation.
    double effective_radius = 1.0;  // Effective radius.
    // Interpolated value.
    auto interpolated_values =
        PUApproximant<Vec2d>::evaluate(domain, values, tree, effective_radius, approx);

    for (auto value : interpolated_values) {
        EXPECT_NEAR(1, value, 1e-15);
    }
}

TEST(Interpolants, PUApproximantKDTree3D) {
    // Domain.
    BoxShape<Vec3d> b(0, 1);
    DomainDiscretization<Vec3d> domain = b.discretizeWithStep(0.1);
    int support_size = 30;
    domain.findSupport(FindClosest(support_size));

    // Values.
    Eigen::VectorXd values(domain.size());
    values.setOnes();

    // Approximation engine.
    RBFFD<Polyharmonic<double>, Vec3d, ScaleToClosest> approx({3}, Monomials<Vec3d>(2));

    // Test points.
    BoxShape<Vec3d> b_test(0, 1);
    DomainDiscretization<Vec3d> domain_test = b_test.discretizeWithStep(0.2);
    auto test_points = domain_test.positions();
    KDTree<Vec3d> tree(test_points);

    // Partition of unity approximation.
    double effective_radius = 1.0;  // Effective radius.
    // Interpolated value.
    auto interpolated_values =
        PUApproximant<Vec3d>::evaluate(domain, values, tree, effective_radius, approx);

    for (auto value : interpolated_values) {
        EXPECT_NEAR(1, value, 1e-15);
    }
}

}  // namespace mm
