#include "medusa/bits/domains/PolyhedronShape.hpp"
#include "gtest/gtest.h"

namespace mm {

TEST(Domains, PolyhedronShapeUsageExample) {
    /// [PolyhedronShape usage example]
    auto poly = PolyhedronShape<Vec3d>::fromOFF("test/testdata/bunny.off");
    if (poly.contains({2.3, 4.5, 0.1})) {
        // do something
    }
    std::cout << poly << std::endl;
    // Boundary types are derived from face colors.
    auto d = poly.discretizeBoundaryWithStep(2.5);
    /// [PolyhedronShape usage example]
    (void) d;
}

TEST(Domains, PolyhedronShapeContains) {
    auto poly = PolyhedronShape<Vec3d>::fromOFF("test/testdata/tetrahedron.off");
    EXPECT_TRUE(poly.contains({0, 0, 0}));
    EXPECT_FALSE(poly.contains({-1, -1, -1}));
    EXPECT_FALSE(poly.contains({-0.1, 0, 0}));
    EXPECT_FALSE(poly.contains({0, -1, 0}));
    EXPECT_FALSE(poly.contains({0, 0, -0.1}));
    EXPECT_TRUE(poly.contains({0.2, 0.2, 0.2}));
}

TEST(Domains, PolyhedronShapeDiscretizeWithType) {
    auto poly = PolyhedronShape<Vec3d>::fromOFF("test/testdata/tetrahedron.off");
    auto d = poly.discretizeBoundaryWithStep(0.1);

    std::vector<int> expected = {
        poly.rgb2type(255, 0, 0),
        poly.rgb2type(0, 255, 0),
        poly.rgb2type(0, 0, 255),
        poly.rgb2type(0, 0, 0),
    };

    for (int type : d.types()) {
        EXPECT_NE(std::find(expected.begin(), expected.end(), type), expected.end());
    }

    auto contains = [&](const Vec3d& p) {
        return std::find(d.positions().begin(), d.positions().end(), p) != d.positions().end();
    };
    EXPECT_TRUE(contains({0.0, 0.0, 0.0}));
    EXPECT_TRUE(contains({1.0, 0.0, 0.0}));
    EXPECT_TRUE(contains({0.0, 1.0, 0.0}));
    EXPECT_TRUE(contains({0.0, 0.0, 1.0}));

    // Size is dependent on the order of triangles and their boundary points.
    EXPECT_EQ(d.size(), 110);
}

TEST(Domains, PolyhedronShapeOnlyTriangular) {
    EXPECT_DEATH(PolyhedronShape<Vec3d>::fromOFF("test/testdata/cube.off"),
                 "Surface must be a triangle mesh");
}

TEST(Domains, PolyhderonToPoltopeShapeCast) {
    auto poly = PolytopeShape<Vec3d>::fromOFF("test/testdata/bunny.off");
    Vec3d bot = {-95.3614, -59.3019, 33.3985};
    Vec3d top = {61.5441, 61.8599, 187.632};
    auto bb = poly.bbox();
    EXPECT_LT((bb.first - bot).norm(), 1e-5);
    EXPECT_LT((bb.second - top).norm(), 1e-5);
}

}  // namespace mm

