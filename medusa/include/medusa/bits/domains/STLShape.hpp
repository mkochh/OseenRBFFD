#ifndef MEDUSA_BITS_DOMAINS_STLSHAPE_HPP_
#define MEDUSA_BITS_DOMAINS_STLSHAPE_HPP_

#include <medusa/Config.hpp>
#include "STLShape_fwd.hpp"
#include "DomainDiscretization_fwd.hpp"
#include "discretization_helpers_advanced.hpp"
#include <medusa/bits/spatial_search/KDTreeMutable.hpp>

#include <map>

/**
 * @file
 * Implementation of STL shapes.
 */

namespace mm {

template <typename vec_t>
STLShape<vec_t>::STLShape(const std::vector<STL::Triangle>& stl) {
    assert_msg(!stl.empty(), "Triangle list should not be empty.");
    int num_triangles = stl.size();
    std::map<vec_t, int> index;
    int num_points = 0;
    faces_.resize(num_triangles);
    normals_.resize(num_triangles);
    bbox_.first = bbox_.second = v(stl[0].p1);
    for (int i = 0; i < num_triangles; ++i) {
        const auto& t = stl[i];

        auto p1 = v(t.p1);
        if (index.count(p1) == 0) index[p1] = num_points++;
        faces_[i][0] = index[p1];
        bbox_.first = bbox_.first.cwiseMin(p1);
        bbox_.second = bbox_.second.cwiseMax(p1);

        auto p2 = v(t.p2);
        if (index.count(p2) == 0) index[p2] = num_points++;
        faces_[i][1] = index[p2];
        bbox_.first = bbox_.first.cwiseMin(p2);
        bbox_.second = bbox_.second.cwiseMax(p2);

        auto p3 = v(t.p3);
        if (index.count(p3) == 0) index[p3] = num_points++;
        faces_[i][2] = index[p3];
        bbox_.first = bbox_.first.cwiseMin(p3);
        bbox_.second = bbox_.second.cwiseMax(p3);

        normals_[i] = v(t.normal).normalized();
    }

    vertices_.resize(num_points);
    for (const auto& p : index) {
        vertices_[p.second] = p.first;
    }
}

template <typename vec_t>
DomainDiscretization<vec_t>
STLShape<vec_t>::discretizeBoundaryWithDensity(
        const std::function<scalar_t(vec_t)>& dx, int type) const {
    if (type == 0) type = -1;
    DomainDiscretization<vec_t> d(*this);

    KDTreeMutable<vec_t> tree;
    int num_faces = faces_.size();
    for (int i = 0; i < num_faces; ++i) {
        auto p1 = vertices_[faces_[i][0]];
        auto p2 = vertices_[faces_[i][1]];
        auto p3 = vertices_[faces_[i][2]];
        auto n = normals_[i];

        Range<vec_t> points = discretization_helpers::discretizeTriangleWithDensity(
                p1, p2, p3, n, dx, false);

        for (const auto& p : points) {
            if (tree.existsPointInSphere(p, dx(p))) continue;
            tree.insert(p);
            d.addBoundaryNode(p, type, n);
        }
    }

    return d;
}

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_STLSHAPE_HPP_
