#include <medusa/bits/domains/GeneralFill_fwd.hpp>
#include <medusa/bits/types/Vec.hpp>
#include <medusa/bits/domains/BallShape.hpp>
#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/domains/DomainDiscretization_fwd.hpp>
#include <medusa/bits/domains/FindClosest.hpp>
#include <medusa/bits/spatial_search/KDGrid.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(DomainEngines, GeneralFill1D) {
    GeneralFill<Vec1d> fill_engine;
    fill_engine.seed(1);
    // create test domain
    BoxShape<Vec1d> b(0.0, 1.0);
    double dx = 0.01;
    DomainDiscretization<Vec1d> domain = b.discretizeBoundaryWithStep(dx);
    // execute fill domain with constant density
    fill_engine(domain, dx);

    // find minimal distance
    domain.findSupport(FindClosest(2));
    std::vector<double> tmp;
    for (int i = 0; i < domain.size(); ++i) {
        tmp.push_back(domain.dr(i));
    }
    // check if minimal distance between two nodes is above target radius * proximityTolerance
    EXPECT_GE(*std::min_element(tmp.begin(), tmp.end()), dx-1e-9);

    domain.findSupport(FindClosest(2));
    for (int i : domain.interior()) {
        double dr = domain.dr(i);
        double err = std::abs(dr - dx) / dx;
        EXPECT_LT(err, 0.31);
    }
}

TEST(DomainEngines, GeneralFill2DConstantFill) {
    GeneralFill<Vec2d> fill_engine;
    // setup
    fill_engine.seed(1);
    // create test domain
    BoxShape<Vec2d> b({0, 0}, {1.2, 1.5});
    BallShape<Vec2d> c({0.4, 0.4}, 0.25);
    DomainDiscretization<Vec2d> domain(b-c);
    // execute fill domain with constant density
    double dx = 0.01;
    fill_engine(domain, dx);
    // find minimal distance
    domain.findSupport(FindClosest(2));
    std::vector<double> tmp;
    for (int i = 0; i < domain.size(); ++i) {
        tmp.push_back(domain.dr(i));
    }
    // check if minimal distance between two nodes is above target radius
    EXPECT_GE(*std::min_element(tmp.begin(), tmp.end()), dx-1e-9);

    domain.findSupport(FindClosest(2));
    for (int i : domain.interior()) {
        double dr = domain.dr(i);
        double err = std::abs(dr - dx) / dx;
        EXPECT_LT(err, 0.2);
    }
}

TEST(DomainEngines, GeneralFill2DFunctionFill) {
    GeneralFill<Vec2d> fill_engine;
    // setup
    fill_engine.seed(1);
    // create test domain
    BoxShape<Vec2d> b({0, 0}, {1.2, 1.5});
    BallShape<Vec2d> c({0.4, 0.4}, 0.25);
    DomainDiscretization<Vec2d> domain(b-c);
    // define target density
    auto dx = [](const Vec2d& p) { return std::pow(p[0]/10, 2) + 0.01; };
    fill_engine(domain, dx);
    // find minimal distance
    domain.findSupport(FindClosest(2));
    std::vector<double> tmp;
    for (int i = 0; i < domain.size(); ++i) {
        tmp.push_back(domain.dr(i));
    }
    // check if minimal distance between two nodes is above target radius
    EXPECT_GE(*std::min_element(tmp.begin(), tmp.end()),  dx(0.0)-1e-9);

    domain.findSupport(FindClosest(2));
    for (int i : domain.interior()) {
        double dr = domain.dr(i);
        double h = dx(domain.pos(i));
        double err = std::abs(dr - h) / h;
        EXPECT_LT(err, 0.25);
    }
}

TEST(DomainEngines, GeneralFill2DFunctionFillBackgroundGrid) {
    GeneralFill<Vec2d> fill_engine;
    // setup
    fill_engine.seed(1);
    // create test domain
    Vec2d bot = {0, 0}, top = {1.2, 1.5};
    BoxShape<Vec2d> b(bot, top);
    BallShape<Vec2d> c({0.4, 0.4}, 0.25);
    DomainDiscretization<Vec2d> domain(b-c);
    // define target density
    auto dx = [](const Vec2d& p) { return std::pow(p[0]/10, 2) + 0.01; };

    /// [Background grid]
    KDGrid<Vec2d> grid(bot, top, 0.005);
    fill_engine(domain, dx, grid);
    /// [Background grid]
    // find minimal distance
    domain.findSupport(FindClosest(2));
    std::vector<double> tmp;
    for (int i = 0; i < domain.size(); ++i) {
        tmp.push_back(domain.dr(i));
    }
    // check if minimal distance between two nodes is above target radius
    EXPECT_GE(*std::min_element(tmp.begin(), tmp.end()),  dx(0.0)-1e-9);

    domain.findSupport(FindClosest(2));
    for (int i : domain.interior()) {
        double dr = domain.dr(i);
        double h = dx(domain.pos(i));
        double err = std::abs(dr - h) / h;
        EXPECT_LT(err, 0.25);
    }
}

TEST(DomainEngines, GeneralFill3DConstantFill) {
    GeneralFill<Vec3d> fill_engine;
    fill_engine.seed(1);
    // create test domain
    BoxShape<Vec3d> b({0, 0, 0}, {1.2, 1.5, 1.1});
    BallShape<Vec3d> c({0, 0, 0}, 0.5);
    DomainDiscretization<Vec3d> domain(b-c);
    // execute fill domain with constant density
    double dx = 0.05;
    fill_engine(domain, dx);

    // find minimal distance - Actual test
    domain.findSupport(FindClosest(2));
    std::vector<double> tmp;
    for (int i = 0; i < domain.size(); ++i) {
        tmp.push_back(domain.dr(i));
    }
    // check if minimal distance between two nodes is above target radius
    EXPECT_GE(*std::min_element(tmp.begin(), tmp.end()), dx-1e-9);

    domain.findSupport(FindClosest(2));
    for (int i : domain.interior()) {
        double dr = domain.dr(i);
        double err = std::abs(dr - dx) / dx;
        EXPECT_LT(err, 0.31);
    }
}

TEST(DomainEngines, GeneralFill3DFunctionFill) {
    /// [PDS usage example]
    GeneralFill<Vec3d> fill_engine;
    fill_engine.seed(1);
    // create test domain
    BoxShape<Vec3d> b({0, 0, 0}, {1.2, 1.5, 1.1});
    BallShape<Vec3d> c({0, 0, 0}, 0.5);
    DomainDiscretization<Vec3d> domain(b-c);
    // define target density
    auto dx = [](const Vec3d& p) {
        return std::pow(p[0]/10 + p[1]/10 + p[2]/10, 2) + 0.025;
    };
    domain.fill(fill_engine, dx);
    /// [PDS usage example]

    // find minimal distance - Actual test
    domain.findSupport(FindClosest(2));
    std::vector<double> tmp;
    for (int i = 0; i < domain.size(); ++i) {
        tmp.push_back(domain.dr(i));
    }
    // check if minimal distance between two nodes is above target radius
    EXPECT_GE(*std::min_element(tmp.begin(), tmp.end()), dx(0.0)-1e-9);

    domain.findSupport(FindClosest(2));
    for (int i : domain.interior()) {
        double dr = domain.dr(i);
        double h = dx(domain.pos(i));
        double err = std::abs(dr - h) / h;
        EXPECT_LT(err, 0.4);
    }
}

TEST(DomainEngines, GeneralFillMaxPoints) {
    GeneralFill<Vec2d> fill_engine;
    int max_points = 1000, num_samples = 15;
    // setup
    fill_engine.seed(1).maxPoints(max_points).numSamples(num_samples);
    // create test domain
    BoxShape<Vec2d> b({0, 0}, {1.2, 1.5});
    BallShape<Vec2d> c({0.4, 0.4}, 0.25);
    DomainDiscretization<Vec2d> domain(b-c);
    // execute fill domain with constant density
    double dx = 0.01;
    fill_engine(domain, dx);

    EXPECT_LE(domain.positions().size(), max_points + num_samples);
    EXPECT_GE(domain.positions().size(), max_points - num_samples);
}

TEST(DomainEngines, GeneralFillAlternativeContains) {
    int seed = 1, max_points = 1000;
    double h = 0.2;

    BallShape<Vec3d> shape({0, 0, 0}, 1);
    GeneralFill<Vec3d> gf;
    gf.seed(seed).maxPoints(max_points);

    auto contains = [&] (Vec3d p) {
        return shape.contains(p) && p(0) > 0;
    };

    KDTreeMutable<Vec3d> tree;
    DomainDiscretization<Vec3d> domain(shape);
    gf(domain, h, tree, contains);

    for (const auto &pt : domain.positions()) {
        EXPECT_GT(pt(0), 0);
    }
}

}  // namespace mm
