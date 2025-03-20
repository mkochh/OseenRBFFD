#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/domains/FindClosest.hpp>
#include "medusa/bits/domains/GrainDropFill.hpp"

#include "gtest/gtest.h"

namespace mm {


TEST(DomainEngines, GrainDrop1d) {
    BoxShape<Vec1d> b(0.0, 1.0);
    double dx = 0.05;
    DomainDiscretization<Vec1d> domain = b.discretizeBoundaryWithStep(dx);
    // execute fill domain with constant density
    GrainDropFill<Vec1d> fill_engine(0.0, 1.0);
    fill_engine.seed(0).initialSpacing(dx/5).maxPoints(3000);
    fill_engine(domain, dx);

    // tests mostly that it compiles and runs, as this invokes a specialization.

    domain.findSupport(FindClosest(2));
    for (int i : domain.interior()) {
        double dr = domain.dr(i);
        double err = std::abs(dr - dx) / dx;
        EXPECT_LT(err, 0.1);
    }
}

TEST(DomainEngines, GrainDrop2d) {
    BoxShape<Vec2d> b(0.0, 1.0);
    auto dx = [](const Vec2d& p) { return std::pow(p[0]/10, 2) + 0.01; };
    DomainDiscretization<Vec2d> domain = b.discretizeBoundaryWithDensity(dx);
    // execute fill domain with constant density
    GrainDropFill<Vec2d> fill_engine(0.0, 1.0);
    fill_engine.seed(0).initialSpacing(0.01/5).maxPoints(15000).excessFactor(10);
    fill_engine(domain, dx);

    domain.findSupport(FindClosest(2));
    for (int i : domain.interior()) {
        double dr = domain.dr(i);
        double h = dx(domain.pos(i));
        double err = std::abs(dr - h) / h;
        EXPECT_LT(err, 0.3);
    }
}

TEST(DomainEngines, GrainDrop3d) {
    /// [GDF usage example]
    BoxShape<Vec3d> b(0.0, 1.0);
    auto dx = [](const Vec3d& p) {
        return std::pow(p[0]/10 + p[1]/10 + p[2]/10, 2) + 0.025;
    };
    DomainDiscretization<Vec3d> domain = b.discretizeBoundaryWithDensity(dx);
    // execute fill domain with constant density
    GrainDropFill<Vec3d> fill_engine(0.0, 1.0);
    fill_engine.seed(0).initialSpacing(0.025/5).maxPoints(50000);
    fill_engine(domain, dx);
    /// [GDF usage example]

    domain.findSupport(FindClosest(2));
    for (int i : domain.interior()) {
        double dr = domain.dr(i);
        double h = dx(domain.pos(i));
        double err = std::abs(dr - h) / h;
        EXPECT_LT(err, 0.4);
    }
}

}  // namespace mm
