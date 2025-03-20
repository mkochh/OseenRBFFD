#include "medusa/bits/domains/PolytopeShape.hpp"
#include "gtest/gtest.h"

namespace mm {

TEST(Domains, PolytopeShapeUsageExample) {
    /// [PolytopeShape usage example]
    // You must include PolyhedronShape.hpp (which depends on CGAL) for the 3D version.
    PolytopeShape<Vec2d> poly({{0.0, -2.0}, {1.0, -1.0}, {2.0, 0.0}, {1.0, 1.0}, {0.0, 2.0},
                               {-1.0, 1.0}, {-2.0, 0.0}, {-1.0, -1.0}});
    if (poly.contains({2.3, 4.5})) {
        // do something
    }
    std::cout << poly << std::endl;
    auto d = poly.discretizeBoundaryWithStep(0.1);
    /// [PolytopeShape usage example]
    (void) d;
}

TEST(Domains, PolytopeShapeCast) {
    PolygonShape<Vec2d> poly({{0.0, -2.0}, {1.0, -1.0}, {2.0, 0.0}, {1.0, 1.0}, {0.0, 2.0},
                              {-1.0, 1.0}, {-2.0, 0.0}, {-1.0, -1.0}});
    // Check that this compiles.
    PolytopeShape<Vec2d> poly2 = poly;
    EXPECT_EQ(poly.points(), poly2.points());
}

}  // namespace mm

