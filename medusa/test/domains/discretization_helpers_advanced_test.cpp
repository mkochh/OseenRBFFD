#include <medusa/bits/domains/discretization_helpers_advanced.hpp>
#include <medusa/bits/types/Vec_fwd.hpp>

#include "gtest/gtest.h"
#include <algorithm>

namespace mm {

TEST(Domains, DiscretizeTriangleWithDensity) {
    /// [discretizeTriangleWithDensity]
    Vec3d p1 = {0, 0, 1}, p2 = {0, 1, 0}, p3 = {1, 0, 0};
    Vec3d n = (p2-p1).cross(p3-p1).normalized();
    auto fn = [](const Vec3d&) { return 0.3; };
    std::vector<Vec3d> pts_with_bnd = discretization_helpers::discretizeTriangleWithDensity(
            p1, p2, p3, n, fn, false);  // include boundary nodes
    std::vector<Vec3d> pts_only_int = discretization_helpers::discretizeTriangleWithDensity(
            p1, p2, p3, n, fn, true);
    /// [discretizeTriangleWithDensity]
    // Plot the points to check visually
    // Interior is a subset of with boundary.
    for (const auto& p : pts_only_int) {
        auto it = std::find(pts_with_bnd.begin(), pts_with_bnd.end(), p);
        EXPECT_NE(it, pts_with_bnd.end());
    }
    // Check the normals of all triplets -- they should all be the same up to a sign change.
    if (n[0] < 0) { n = -n; }
    int s = pts_with_bnd.size();
    for (int i = 0; i < s; ++i) {
        for (int j = i+1; j < s; ++j) {
            auto q1 = pts_with_bnd[i], q2 = pts_with_bnd[j];
            for (int k = j+1; k < s; ++k) {
                auto q3 = pts_with_bnd[k];
                Vec3d local_normal = (q1-q2).cross(q1-q3);
                double len = local_normal.norm();
                if (len < 1e-14) continue;  // skip colinear points
                local_normal /= len;
                if (local_normal[0] < 0) { local_normal = -local_normal; }
                EXPECT_LT((local_normal-n).norm(), 2e-14);
            }
            // Check the distance. It can get closer than dx near the corners.
            EXPECT_GE((q2-q1).norm(), 0.2);
        }
    }
}

}  // namespace mm
