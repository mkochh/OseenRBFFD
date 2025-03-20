#ifndef MEDUSA_BITS_DOMAINS_BALLSHAPE_HPP_
#define MEDUSA_BITS_DOMAINS_BALLSHAPE_HPP_

#include "BallShape_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>
#include "DomainDiscretization.hpp"
#include "GeneralFill.hpp"
#include "discretization_helpers.hpp"
#include <medusa/bits/utils/numutils.hpp>

/**
 * @file
 * Implementation of class for ball shaped domains.
 */

namespace mm {

template <typename vec_t>
DomainDiscretization<vec_t> BallShape<vec_t>::discretizeBoundaryWithStep(scalar_t step,
                                                                         int type) const {
    if (type == 0) type = -1;
    DomainDiscretization<vec_t> domain(*this);
    int num = iceil(2*Pi<scalar_t>::value*radius_/step);
    auto points = discretization_helpers::SphereDiscretization<scalar_t, dim>::construct(
            radius_, num);
    for (const auto& p : points) {
        domain.addBoundaryNode(center_+p, type, p.normalized());
    }
    return domain;
}

template <typename vec_t>
DomainDiscretization<vec_t> BallShape<vec_t>::discretizeWithStep(
        scalar_t step, int internal_type, int boundary_type) const {
    auto d = discretizeBoundaryWithStep(step, boundary_type);
    if (internal_type == 0) internal_type = 1;
    if (dim == 1) {
        int num_of_points = iceil(2*radius_/step) - 1;
        for (const vec_t& v : linspace(bbox().first, bbox().second, {num_of_points}, false)) {
            d.addInternalNode(v, internal_type);
        }
    } else if (dim == 2) {
        // Use concentric circles with step dx.
        int num_of_circles = iceil(radius_ / step);  // not counting outer circle
        scalar_t dr = radius_ / num_of_circles;
        d.addInternalNode(center_, internal_type);
        for (int i = 1; i < num_of_circles; ++i) {
            scalar_t r = i*dr;
            int num_of_points = iceil(Pi<scalar_t>::value / std::asin(step/2.0/r));
            scalar_t dfi = 2*Pi<scalar_t>::value / num_of_points;
            for (int j = 0; j < num_of_points; ++j) {
                d.addInternalNode(center_ + r * vec_t({std::cos(j*dfi), std::sin(j*dfi)}),
                                  internal_type);
            }
        }
    } else if (dim >= 3) {  // Call default fill engine (could do concentric n-dim spheres).
        assert_msg(false, "This domain does not support filling with density, "
                          "use a GeneralSurfaceFill engine instead.");
    }
    return d;
}

template <typename vec_t>
std::ostream& BallShape<vec_t>::print(std::ostream& os) const {
    return os << "BallShape(" << center_.transpose() << ", " << radius_ << ")";
}

template <typename vec_t>
DomainDiscretization<vec_t>
BallShape<vec_t>::discretizeBoundaryWithDensity(const std::function<scalar_t(vec_t)>& dr,
                                                int type) const {
    if (type == 0) type = -1;
    DomainDiscretization<vec_t> domain(*this);
    if (dim == 1) {  // always only two points
        domain.addBoundaryNode({center_[0] - radius_}, type, vec_t(-1.0));
        domain.addBoundaryNode({center_[0] + radius_}, type, vec_t(1.0));
    } else if (dim == 2) {
        scalar_t cur_alpha = 0.0;
        scalar_t max_alpha = 2*Pi<scalar_t>::value;
        Range<scalar_t> alphas = {cur_alpha};
        while (cur_alpha < max_alpha) {
            scalar_t r = dr(center_ + radius_ * vec_t({std::cos(cur_alpha), std::sin(cur_alpha)}));
            cur_alpha += 2*std::asin(r/2.0/radius_);
            alphas.push_back(cur_alpha);
        }

        scalar_t shrink = 2*Pi<scalar_t>::value / cur_alpha;
        for (int i = 0; i < alphas.size()-1; ++i) {
            scalar_t a = alphas[i] * shrink;
            vec_t normal = {std::cos(a), std::sin(a)};
            domain.addBoundaryNode(center_ + radius_ * normal, type, normal);
        }
    } else if (dim == 3) {
        // define modified density
        auto dr2upper = [&](const Vec<scalar_t, 2>& v) {
            double d = (1 + v.squaredNorm());
            vec_t p = {2*v[0]/d, 2*v[1]/d, (1-v.squaredNorm())/d};
            return d/2 * dr(center_ + radius_ * p) / radius_;
        };

        BallShape<Vec<scalar_t, 2>> circ(0.0, 1.0);
        GeneralFill<Vec<scalar_t, 2>> fill;
        fill.seed(0).proximityTolerance(0.99);
        auto domain2d = circ.discretizeWithDensity(dr2upper, fill);

        for (int i : domain2d.interior()) {
            auto v = domain2d.pos(i);
            double d = (1 + v.squaredNorm());
            vec_t p = {2*v[0]/d, 2*v[1]/d, (1-v.squaredNorm())/d};
            domain.addBoundaryNode(center_ + radius_*p, type, p.normalized());
        }

        for (int i : domain2d.boundary()) {
            auto v = domain2d.pos(i);
            vec_t p = {v[0], v[1], 0};
            domain.addBoundaryNode(center_ + radius_*p, type, p.normalized());
        }

        // lower hemisphere
        auto dr2lower = [&](const Vec<scalar_t, 2>& v) {
            double d = (1 + v.squaredNorm());
            vec_t p = {2*v[0]/d, 2*v[1]/d, -(1-v.squaredNorm())/d};
            return d/2 * dr(center_ + radius_ * p) / radius_;
        };
        domain2d = circ.discretizeWithDensity(dr2lower, fill);
        for (int i : domain2d.interior()) {
            auto v = domain2d.pos(i);
            double d = (1 + v.squaredNorm());
            vec_t p = {2*v[0]/d, 2*v[1]/d, -(1-v.squaredNorm())/d};
            domain.addBoundaryNode(center_ + radius_*p, type, p.normalized());
        }
    } else {
        assert_msg(false, "This domain does not support filling with density, "
                          "use a GeneralSurfaceFill engine instead.");
    }
    return domain;
}

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_BALLSHAPE_HPP_
