#ifndef MEDUSA_BITS_DOMAINS_DISCRETIZATION_HELPERS_ADVANCED_HPP_
#define MEDUSA_BITS_DOMAINS_DISCRETIZATION_HELPERS_ADVANCED_HPP_

/**
 * @file
 * Contains helpers for discretizing various geometric objects that need more advanced
 * parts of Medusa to work.
 *
 * @example test/domains/discretization_helpers_advanced_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/domains/PolygonShape.hpp>
#include <medusa/bits/domains/GeneralFill.hpp>
#include <medusa/bits/types/Range.hpp>

#include <Eigen/Geometry>

#include <vector>

namespace mm {
namespace discretization_helpers {

/**
 * Discretize a triangle in 3D space with given density.
 *
 * The orientation of the triangle is not important. The triangle is discretized by first
 * discretizing the boundary and filling the interior with a GeneralFill.
 *
 * @tparam vec_t Must be a three-dimensional vector type.
 * @param p1 First point.
 * @param p2 Second point.
 * @param p3 Third point.
 * @param normal Unit normal vector of the plane defined by the triangle.
 * @param h Spacing function.
 * @param only_interior If true, only the points in the interior are returned, otherwise, boundary
 * discretization points are returned as well.
 * @return A list of discretization points in the triangle.
 *
 * @ingroup domains
 * @snippet domains/discretization_helpers_advanced_test.cpp discretizeTriangleWithDensity
 * @sa GeneralFill, PolygonShape, STLShape
 */
template <typename vec_t, typename func_t>
Range<vec_t> discretizeTriangleWithDensity(const vec_t& p1, const vec_t& p2, const vec_t& p3,
                                           const vec_t& normal, const func_t& h,
                                           bool only_interior = true) {
    static_assert(vec_t::dim == 3, "This method is intended for 3D vectors only.");
    using scalar_t = typename vec_t::scalar_t;
    using V2 = mm::Vec<scalar_t, 2>;
    static constexpr int dim = vec_t::dim;

    // Rotate to flat
    vec_t zunit = {0, 0, 1};
    scalar_t c = normal.dot(zunit);
    if (c < 0) {
        zunit = -zunit;
        c = -c;
    }
    vec_t vp = normal.cross(zunit);

    Eigen::Matrix<scalar_t, dim, dim> R;
    R.setIdentity();
    Eigen::Matrix<scalar_t, dim, dim> VP;
    VP << 0, -vp[2], vp[1], vp[2], 0, -vp[0], -vp[1], vp[0], 0;
    R += VP + 1 / (1 + c) * VP * VP;

    vec_t np1 = R * p1;
    vec_t np2 = R * p2;
    vec_t np3 = R * p3;
    // All z coordinates are the same, but still take the average.
    scalar_t z = 1.0 / 3.0 * (np1[2] + np2[2] + np3[2]);


    // Fill the flat triangle.
    std::vector<V2> pts = {np1.template head<2>(), np2.template head<2>(),
                           np3.template head<2>()};
    mm::PolygonShape<V2> triangle(pts);
    auto new_h = [&](const V2& p) { return h(R.transpose() * vec_t(p[0], p[1], z)); };
    auto domain = triangle.discretizeBoundaryWithDensity(new_h);
    mm::GeneralFill<V2> fill;
    fill.seed(0);
    fill(domain, new_h);

    auto idxs = only_interior ? domain.interior() : domain.all();
    std::vector<vec_t> result;
    result.reserve(idxs.size());
    for (int i : idxs) {
        auto ip = domain.pos(i);
        result.emplace_back(R.transpose() * vec_t(ip[0], ip[1], z));
    }
    return result;
}

}  // namespace discretization_helpers
}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_DISCRETIZATION_HELPERS_ADVANCED_HPP_
