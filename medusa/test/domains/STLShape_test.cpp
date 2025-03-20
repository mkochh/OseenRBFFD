#include "medusa/bits/domains/STLShape.hpp"
#include "medusa/bits/domains/FindClosest.hpp"
#include "medusa/bits/domains/DomainDiscretization_fwd.hpp"

#include "gtest/gtest.h"

namespace mm {

TEST(Domains, STLShape) {
    /// [STL shape usage example]
    auto data = STL::read("test/testdata/tetrahedron.stl");
    STLShape<Vec3d> shape(data);
    auto dx = [](const Vec3d& p) { return 0.005 + (p-Vec3d(0.5, 0.5, 0.5)).norm()/50; };
    auto d = shape.discretizeBoundaryWithDensity(dx);
    /// [STL shape usage example]

    EXPECT_GT(d.size(), 5000);

    d.findSupport(FindClosest(2));

    for (int i = 0; i < d.size(); ++i) {
        ASSERT_NEAR(d.dr(i), dx(d.pos(i)), 2e-2);
    }

    auto bbox = shape.bbox();
    EXPECT_EQ(Vec3d(0.5, 0.5, 0.5), bbox.first);
    EXPECT_EQ(Vec3d(1.5, 1.5, 1.5), bbox.second);
}

}  // namespace mm
