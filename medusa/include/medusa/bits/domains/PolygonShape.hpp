#ifndef MEDUSA_BITS_DOMAINS_POLYGONSHAPE_HPP_
#define MEDUSA_BITS_DOMAINS_POLYGONSHAPE_HPP_

#include "PolygonShape_fwd.hpp"
#include "discretization_helpers.hpp"
#include "DomainDiscretization.hpp"
#include <medusa/bits/utils/numutils.hpp>

/**
 * @file
 * Implementation of PolygonShape class.
 */

namespace mm {

template <typename vec_t>
PolygonShape<vec_t>::PolygonShape(const std::vector<vec_t>& points) : points_(points) {
    scalar_t area = 0;
    int n = points_.size();
    assert_msg(n >= 3, "At least three points are needed to form a polygon, got %d.", n);
    for (int i = 0; i < n; ++i) {
        int j = (i == n-1) ? 0 : (i+1);  // next point
        area += (points_[i][0] - points_[j][0]) * (points_[i][1] + points_[j][1]);
    }
    assert_msg(std::abs(area) > 1e-10, "Given polygon has no area.");
    if (area < 0) {
        std::reverse(points_.begin(), points_.end());
    }
    extendByMargin();
}

template <typename vec_t>
void PolygonShape<vec_t>::setMargin(scalar_t margin) {
    base_t::setMargin(margin);
    extendByMargin();
}

template <typename vec_t>
void PolygonShape<vec_t>::extendByMargin() {
    int n = points_.size();
    points_with_margin_ = points_;
    for (int i = 0; i < n; ++i) {
        int j = (i == n - 1) ? 0 : (i + 1);  // next
        int k = (i == 0) ? n - 1 : (i - 1);  // prev
        vec_t next_edge = (points_[j] - points_[i]).normalized();
        vec_t prev_edge = (points_[i] - points_[k]).normalized();
        vec_t normal = {next_edge[1] + prev_edge[1], -next_edge[0] - prev_edge[0]};
        points_with_margin_[i] += margin_*normal.normalized();
    }
    points_with_margin_.push_back(points_with_margin_[0]);
}

template <typename vec_t>
bool PolygonShape<vec_t>::contains(const vec_t& point) const {
    // Winding number test loosely based on http://geomalgorithms.com/a03-_inclusion.html.
    int n = points_.size();  // loop through all edges of the polygon
    int wn = 0;  // the winding number counter
    auto& pts = points_with_margin_;
    for (int i = 0; i < n; i++) {  // edge from points[i] to points[i+1]
        if (pts[i][1] <= point[1]) {  // Start y <= point[1].
            // An upward crossing and point left of edge.
            if (pts[i+1][1] > point[1] && isLeft(pts[i], pts[i+1], point) > 0) {
                ++wn;  // Have a valid up intersect.
            }
        } else {  // Start y > point[1] (no test needed).
            // A downward crossing and point right of edge.
            if (pts[i+1][1] <= point[1] && isLeft(pts[i], pts[i+1], point) < 0) {
                --wn;  // Have a valid down intersect.
            }
        }
    }
    return wn != 0;
}

template <typename vec_t>
std::pair<vec_t, vec_t> PolygonShape<vec_t>::bbox() const {
    scalar_t maxx = -INF, maxy = -INF, minx = INF, miny = INF;
    int n = points_.size();
    for (int i = 0; i < n; ++i) {
        if (points_[i][0] > maxx) maxx = points_[i][0];
        if (points_[i][0] < minx) minx = points_[i][0];
        if (points_[i][1] > maxy) maxy = points_[i][1];
        if (points_[i][1] < miny) miny = points_[i][1];
    }
    return {{minx, miny}, {maxx, maxy}};
}

template <typename vec_t>
DomainDiscretization<vec_t> PolygonShape<vec_t>::discretizeBoundaryWithStep(
        scalar_t step, int type) const {
    DomainDiscretization<vec_t> d(*this);
    int n = points_.size();
    for (int i = 0; i < n; ++i) {
        int j = (i == n-1) ? 0 : (i+1);  // next
        int k = (i == 0) ? n-1 : (i-1);  // prev
        vec_t edge = points_[j] - points_[i];
        scalar_t len = edge.norm();

        // Add the vertex point first.
        vec_t next_edge_n = edge / len;
        vec_t prev_edge_n = (points_[i] - points_[k]).normalized();
        vec_t p_normal = {next_edge_n[1] + prev_edge_n[1], -next_edge_n[0] - prev_edge_n[0]};
        d.addBoundaryNode(points_[i], (type == 0) ? -i-1 : type, p_normal.normalized());

        // And then the interior points on the edge.
        int count = iceil(len/step);
        vec_t e_normal = {edge[1] / len, -edge[0] / len};
        for (int c = 1; c < count; ++c) {  // only internal points
            d.addBoundaryNode(points_[i] + c / scalar_t(count) * edge,
                              (type == 0) ? -i-1 : type, e_normal);
        }
    }
    return d;
}

template <typename vec_t>
DomainDiscretization<vec_t>
PolygonShape<vec_t>::discretizeBoundaryWithDensity(
        const std::function<scalar_t(vec_t)>& dr, int type) const {
    DomainDiscretization<vec_t> d(*this);
    int n = points_.size();
    for (int i = 0; i < n; ++i) {
        int j = (i == n-1) ? 0 : (i+1);  // next
        int k = (i == 0) ? n-1 : (i-1);  // prev
        vec_t edge = points_[j] - points_[i];
        scalar_t len = edge.norm();

        // Add the vertex point first.
        vec_t next_edge_n = edge / len;
        vec_t prev_edge_n = (points_[i] - points_[k]).normalized();
        vec_t p_normal = {next_edge_n[1] + prev_edge_n[1], -next_edge_n[0] - prev_edge_n[0]};
        d.addBoundaryNode(points_[i], (type == 0) ? -i-1 : type, p_normal.normalized());

        auto points = discretization_helpers::discretizeLineWithDensity(points_[i], points_[j], dr);
        vec_t e_normal = {points_[j][1] - points_[i][1], -points_[j][0] + points_[i][0]};
        e_normal.normalize();
        for (const auto& p : points) {  // only internal points
            d.addBoundaryNode(p, (type == 0) ? -i-1 : type, e_normal);
        }
    }
    return d;
}

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_POLYGONSHAPE_HPP_
