#include <medusa/bits/domains/BasicRelax.hpp>
#include <medusa/bits/types/Vec.hpp>
#include <medusa/bits/domains/BallShape.hpp>
#include <medusa/bits/domains/FindClosest.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(DomainEngines, RelaxWithConstDistributionTestBoundaryProjection0) {
    double r = 0.25;
    BallShape<Vec2d> sh({0.5, 0.5}, r);
    double step = 1.0 / 55;
    DomainDiscretization<Vec2d> c = sh.discretizeWithStep(step);

    // no of boundary nodes before relax
    int N_1 = c.boundary().size();
    double min_spacing = 0.5 * step;
    BasicRelax relax;
    relax.iterations(20).initialHeat(5).numNeighbours(4).projectionType(BasicRelax::DO_NOT_PROJECT);
    relax(c, min_spacing);
    // no of boundary nodes after relax
    int N_2 = c.boundary().size();
    // actual tests
    c.findSupport(FindClosest(2));
    std::vector<double> tmp;
    for (int i = 0; i < c.size(); ++i) tmp.push_back((c.pos(i) - c.pos(c.support(i, 1))).norm());
    // check for overlapping nodes
    EXPECT_GT(*std::min_element(std::begin(tmp), std::end(tmp)), min_spacing / 10);
    // check if nodes are added on boundary
    EXPECT_EQ(N_2, N_1);
}

TEST(DomainEngines, RelaxWithConstDistributionTestBoundaryProjection1) {
    double r = 0.25;
    BallShape<Vec2d> sh({0.5, 0.5}, r);
    double step = 1.0 / 55.0;
    DomainDiscretization<Vec2d> c = sh.discretizeWithStep(step);
    Range<int> to_remove;
    int j = 0;
    for (int i : c.boundary()) {
        if (j++ % 3 != 0) to_remove.push_back(i);
    }
    c.removeNodes(to_remove);  // keep every fifth node
    // no of boundary nodes before relax
    int N_1 = c.boundary().size();
    double min_spacing = 0.5 * step;
    BasicRelax relax;
    relax.iterations(50).initialHeat(1).numNeighbours(4).finalHeat(0.1)
            .projectionType(BasicRelax::PROJECT_IN_DIRECTION).boundaryProjectionThreshold(0.35);
    relax(c, min_spacing);
    // no of boundary nodes after relax
    int N_2 = c.boundary().size();
    // actual tests
    c.findSupport(FindClosest(2));
    std::vector<double> tmp;
    for (int i = 0; i < c.size(); ++i) tmp.push_back((c.pos(i) - c.pos(c.support(i, 1))).norm());
    // check for overlapping nodes
    EXPECT_GT(*std::min_element(std::begin(tmp), std::end(tmp)), min_spacing / 10);
    // check if nodes are added on boundary
    EXPECT_GT(N_2, N_1);
}

TEST(DomainEngines, RelaxWithConstDistributionTestBoundaryProjection2) {
    double r = 0.25;
    BallShape<Vec2d> sh({0.5, 0.5}, r);
    double step = 1.0 / 55.0;
    DomainDiscretization<Vec2d> c = sh.discretizeWithStep(step);
    Range<int> to_remove;
    int j = 0;
    for (int i : c.boundary()) {
        if (j++ % 3 != 0) to_remove.push_back(i);
    }
    c.removeNodes(to_remove);  // keep every fifth node
    // no of boundary nodes before relax
    int N_1 = c.boundary().size();
    double min_spacing = 0.5 * step;
    BasicRelax relax;
    relax.iterations(50).initialHeat(1).numNeighbours(5)
            .projectionType(BasicRelax::PROJECT_BETWEEN_CLOSEST);
    relax(c, min_spacing);
    // no of boundary nodes after relax
    int N_2 = c.boundary().size();
    // actual tests
    c.findSupport(FindClosest(2));
    std::vector<double> tmp;
    for (int i = 0; i < c.size(); ++i) tmp.push_back((c.pos(i) - c.pos(c.support(i, 1))).norm());
    // check for overlapping nodes
    EXPECT_GT(*std::min_element(std::begin(tmp), std::end(tmp)), min_spacing / 10);
    // check if nodes are added on boundary
    EXPECT_GT(N_2, N_1);
}

TEST(DomainEngines, RelaxWithConstDistributionTestBoundaryProjectiontreshold) {
    double r = 0.25;
    BallShape<Vec2d> sh({0.5, 0.5}, r);
    double step = 1.0 / 55.0;
    DomainDiscretization<Vec2d> c = sh.discretizeWithStep(step);
    // no of boundary nodes before relax
    int N_1 = c.boundary().size();
    double min_spacing = 0.5 * step;
    BasicRelax relax;
    relax.iterations(20).initialHeat(1).numNeighbours(4).boundaryProjectionThreshold(1.5)
            .projectionType(BasicRelax::PROJECT_IN_DIRECTION);
    relax(c, min_spacing);

    // no of boundary nodes after relax
    int N_2 = c.boundary().size();
    // actual tests
    c.findSupport(FindClosest(2));
    std::vector<double> tmp;
    for (int i = 0; i < c.size(); ++i) tmp.push_back((c.pos(i) - c.pos(c.support(i, 1))).norm());
    // check for overlapping nodes
    EXPECT_GT(*std::min_element(std::begin(tmp), std::end(tmp)), min_spacing / 10);
    // check if nodes are added on boundary
    EXPECT_EQ(N_2, N_1);
}

TEST(DomainEngines, RelaxWithVariableDistribution) {
    double r = 0.25;
    BallShape<Vec2d> sh({0.5, 0.5}, r);

    auto fill_density = [](const Vec2d& p) {
        return (0.005 + (p[0] - 0.5) * (p[0] - 0.5) / 2 + (p[1] - 0.5) * (p[1] - 0.5) / 2);
    };

    /// [BasicRelax usage example]
    DomainDiscretization<Vec2d> domain = sh.discretizeBoundaryWithStep(fill_density({r, 0.0}));
    GeneralFill<Vec2d> fill_engine; fill_engine.seed(1);

    BasicRelax relax_engine;
    relax_engine.iterations(50).initialHeat(1).finalHeat(0.05).boundaryProjectionThreshold(0.75)
            .projectionType(BasicRelax::PROJECT_IN_DIRECTION).numNeighbours(5);

    domain.fill(fill_engine, fill_density);
    domain.relax(relax_engine, fill_density);
    /// [BasicRelax usage example]

    domain.findSupport(FindClosest(2));
    std::vector<double> tmp;
    for (int i = 0; i < domain.size(); ++i) tmp.push_back(domain.dr(i));
    // check for overlapping nodes
    EXPECT_GT(*std::min_element(std::begin(tmp), std::end(tmp)), 0.001);
}

TEST(DomainEngines, RelaxDeathTest) {
    double r = 0.25;
    BallShape<Vec2d> sh({0.5, 0.5}, r);

    BasicRelax relax;

    relax.iterations(50).initialHeat(20).finalHeat(0.05).boundaryProjectionThreshold(0.75)
            .projectionType(BasicRelax::PROJECT_IN_DIRECTION).numNeighbours(5);
    auto fill_density = [](Vec2d p) -> double {
        return (0.05 + (p[0] - 0.5) * (p[0] - 0.5) / 2 + (p[1] - 0.5) * (p[1] - 0.5) / 2);
    };

    double step = 1.0 / 15.0;
    DomainDiscretization<Vec2d> c = sh.discretizeWithStep(step);

    // check for overlapping nodes
    EXPECT_DEATH(relax(c, fill_density), "No nodes in relax pool anymore, perhaps use lower heat");
}

TEST(DomainEngines, RelaxWithoutDistribution) {
    double r = 0.25;
    BallShape<Vec2d> sh({0.5, 0.5}, r);

    BasicRelax relax;
    relax.iterations(100).initialHeat(1).finalHeat(0).numNeighbours(3)
            .projectionType(BasicRelax::PROJECT_IN_DIRECTION).boundaryProjectionThreshold(0.55);
    auto fill_density = [](Vec2d p) -> double {
        return (0.005 + (p[0] - 0.5) * (p[0] - 0.5) / 2 + (p[1] - 0.5) * (p[1] - 0.5) / 2);
    };

    DomainDiscretization<Vec2d> c = sh.discretizeBoundaryWithStep(fill_density(Vec2d({r, 0.0})));
    GeneralFill<Vec2d> fill_engine;
    fill_engine.seed(1);
    fill_engine(c, fill_density);
    relax(c);

    c.findSupport(FindClosest(2));
    std::vector<double> tmp;
    for (int i = 0; i < c.size(); ++i) tmp.push_back((c.pos(i) - c.pos(c.support(i, 1))).norm());
    // check for overlapping nodes
    EXPECT_GT(*std::min_element(std::begin(tmp), std::end(tmp)), 0.001);
}

TEST(DomainEngines, DISABLED_RelaxWikiImage1) {
    double r = 0.25;
    BallShape<Vec2d> sh({0.5, 0.5}, r);

    auto fill_density = [](Vec2d p) -> double {
        return (0.0025 + (p[0] - 0.5) * (p[0] - 0.5) / 2 + (p[1] - 0.5) * (p[1] - 0.5) / 4);
    };
    prn(fill_density(Vec2d({0.25, 0})))
    DomainDiscretization<Vec2d> c = sh.discretizeBoundaryWithStep(
            fill_density(Vec2d({0.25, 0})) / 2);
    GeneralFill<Vec2d> fill_engine;
    fill_engine(c, fill_density);

    BasicRelax relax;
    relax.iterations(200).initialHeat(2).numNeighbours(3)
            .projectionType(BasicRelax::PROJECT_BETWEEN_CLOSEST);
    relax.boundaryProjectionThreshold(0.45).finalHeat(0.5);
    relax(c, fill_density);


    c.findSupport(FindClosest(2));
    std::vector<double> tmp;
    for (int i = 0; i < c.size(); ++i) tmp.push_back((c.pos(i) - c.pos(c.support(i, 1))).norm());
    // check for overlapping nodes
    EXPECT_GT(*std::min_element(std::begin(tmp), std::end(tmp)), 0.0005);
}

TEST(DomainEngines, DISABLED_RelaxWikiImage2) {
    double r = 0.25;
    BallShape<Vec2d> sh({0.5, 0.5}, r);

    auto fill_density = [](Vec2d p) -> double {
        return (0.0025 + (p[0] - 0.5) * (p[0] - 0.5) / 2 + (p[1] - 0.5) * (p[1] - 0.5) / 4);
    };
    prn(fill_density(Vec2d({0.25, 0})))

    DomainDiscretization<Vec2d> c = sh.discretizeBoundaryWithStep(
            fill_density(Vec2d({0.25, 0})) / 2);
    GeneralFill<Vec2d> fill_engine;
    fill_engine(c, fill_density);

    BasicRelax relax;
    relax.iterations(200).initialHeat(2).numNeighbours(3)
            .projectionType(BasicRelax::PROJECT_BETWEEN_CLOSEST);
    relax.boundaryProjectionThreshold(0.45).finalHeat(0.5);
    relax(c);


    c.findSupport(FindClosest(2));
    std::vector<double> tmp;
    for (int i = 0; i < c.size(); ++i) tmp.push_back((c.pos(i) - c.pos(c.support(i, 1))).norm());
    // check for overlapping nodes
    EXPECT_GT(*std::min_element(std::begin(tmp), std::end(tmp)), 0.0005);
}

}  // namespace mm
