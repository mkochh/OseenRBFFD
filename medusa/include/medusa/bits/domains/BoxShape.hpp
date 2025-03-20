#ifndef MEDUSA_BITS_DOMAINS_BOXSHAPE_HPP_
#define MEDUSA_BITS_DOMAINS_BOXSHAPE_HPP_

/**
 * @file
 * Implementation of class for box shaped domains.
 */

#include "BoxShape_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>
#include "DomainDiscretization.hpp"
#include "discretization_helpers.hpp"
#include "DomainShape.hpp"
#include "GeneralFill.hpp"
#include <medusa/bits/utils/numutils.hpp>
#include <medusa/bits/types/Vec.hpp>

namespace mm {

template <typename vec_t>
BoxShape<vec_t>::BoxShape(vec_t beg, vec_t end) : beg_(beg), end_(end) {
    for (int i = 0; i < dim; ++i) {
        assert_msg(std::abs(beg_[i] - end_[i]) >= margin_,
                   "Box has no volume, sides in dimension %d have same coordinate %g.",
                   i, beg_[i]);
        if (beg_[i] > end_[i]) std::swap(beg_[i], end_[i]);
    }
}

template <typename vec_t>
bool BoxShape<vec_t>::contains(const vec_t& point) const {
    for (int i = 0; i < dim; ++i) {
        if (!(beg_[i] - margin_ <= point[i] && point[i] <= end_[i] + margin_)) {
            return false;
        }
    }
    return true;
}

template <typename vec_t>
DomainDiscretization<vec_t> BoxShape<vec_t>::discretizeBoundaryWithStep(scalar_t step,
                                                                        int type) const {
    Vec<int, dim> counts;
    for (int i = 0; i < dim; ++i) counts[i] = iceil((end_[i] - beg_[i]) / step) + 1;
    return discretizeBoundary(counts, type);
}

template <typename vec_t>
DomainDiscretization<vec_t> BoxShape<vec_t>::discretizeBoundary(const Vec<int, dim>& counts,
                                                                int type) const {
    for (int x : counts) {
        assert_msg(x >= 2, "All counts must be greater than 2, got counts = %s", counts);
    }
    DomainDiscretization<vec_t> domain(*this);
    vec_t normal, point;
    if (dim == 1) {
        domain.addBoundaryNode(beg_, type == 0 ? -1 : type, vec_t(-1.0));
        domain.addBoundaryNode(end_, type == 0 ? -2 : type, vec_t(1.0));
    } else if (dim == 2) {
        const int TOP = (type == 0) ? -4 : type;
        const int RIGHT = (type == 0) ? -2 : type;
        const int BOTTOM = (type == 0) ? -3 : type;
        const int LEFT = (type == 0) ? -1 : type;

        scalar_t step = (end_[0] - beg_[0]) / (counts[0] - 1);
        for (int i = 1; i < counts[0]-1; ++i) {
            point[0] = beg_[0] + i*step;
            point[1] = beg_[1];
            normal[0] = 0; normal[1] = -1;
            domain.addBoundaryNode(point, BOTTOM, normal);
            point[1] = end_[1];
            normal[0] = 0; normal[1] = 1;
            domain.addBoundaryNode(point, TOP, normal);
        }
        step = (end_[1] - beg_[1]) / (counts[1] - 1);
        for (int i = 1; i < counts[1]-1; ++i) {
            point[1] = beg_[1] + i*step;
            point[0] = beg_[0];
            normal[0] = -1; normal[1] = 0;
            domain.addBoundaryNode(point, LEFT, normal);
            point[0] = end_[0];
            normal[0] = 1; normal[1] = 0;
            domain.addBoundaryNode(point, RIGHT, normal);
        }
        domain.addBoundaryNode(beg_, LEFT, vec_t({-1, -1}).normalized());
        domain.addBoundaryNode(end_, RIGHT, vec_t({1, 1}).normalized());
        domain.addBoundaryNode(vec_t({beg_[0], end_[1]}), TOP, vec_t({-1, 1}).normalized());
        domain.addBoundaryNode(vec_t({end_[0], beg_[1]}), BOTTOM, vec_t({1, -1}).normalized());
    } else if (dim == 3) {
        int LEFT = -1, RIGHT = -2, FRONT = -3, BACK = -4, BOTTOM = -5, TOP = -6;
        if (type != 0) LEFT = RIGHT = FRONT = BACK = BOTTOM = TOP = type;
        scalar_t step = (end_[0] - beg_[0]) / (counts[0] - 1);  // x
        for (int i = 1; i < counts[0]-1; ++i) {
            point[0] = beg_[0] + i*step;
            point[1] = beg_[1];
            point[2] = beg_[2];
            normal = vec_t({0, -1, -1}).normalized();
            domain.addBoundaryNode(point, BOTTOM, normal);
            point[1] = beg_[1];
            point[2] = end_[2];
            normal = vec_t({0, -1, 1}).normalized();
            domain.addBoundaryNode(point, FRONT, normal);
            point[1] = end_[1];
            point[2] = beg_[2];
            normal = vec_t({0, 1, -1}).normalized();
            domain.addBoundaryNode(point, BACK, normal);
            point[1] = end_[1];
            point[2] = end_[2];
            normal = vec_t({0, 1, 1}).normalized();
            domain.addBoundaryNode(point, TOP, normal);

            scalar_t inner_step = (end_[1] - beg_[1]) / (counts[1] - 1);  // y
            for (int j = 1; j < counts[1]-1; ++j) {
                point[1] = beg_[1] + j*inner_step;
                point[2] = beg_[2];
                normal = vec_t({0, 0, -1});
                domain.addBoundaryNode(point, BOTTOM, normal);
                point[2] = end_[2];
                normal = vec_t({0, 0, 1});
                domain.addBoundaryNode(point, TOP, normal);
            }
            inner_step = (end_[2] - beg_[2]) / (counts[2] - 1);  // z
            for (int j = 1; j < counts[2]-1; ++j) {
                point[2] = beg_[2] + j*inner_step;
                point[1] = beg_[1];
                normal = vec_t({0, -1, 0});
                domain.addBoundaryNode(point, FRONT, normal);
                point[1] = end_[1];
                normal = vec_t({0, 1, 0});
                domain.addBoundaryNode(point, BACK, normal);
            }
        }
        step = (end_[1] - beg_[1]) / (counts[1] - 1);  // y
        for (int i = 1; i < counts[1]-1; ++i) {
            point[1] = beg_[1] + i*step;
            point[0] = beg_[0];
            point[2] = beg_[2];
            normal = vec_t({-1, 0, -1}).normalized();
            domain.addBoundaryNode(point, LEFT, normal);
            point[0] = beg_[0];
            point[2] = end_[2];
            normal = vec_t({-1, 0, 1}).normalized();
            domain.addBoundaryNode(point, TOP, normal);
            point[0] = end_[0];
            point[2] = beg_[2];
            normal = vec_t({1, 0, -1}).normalized();
            domain.addBoundaryNode(point, BOTTOM, normal);
            point[0] = end_[0];
            point[2] = end_[2];
            normal = vec_t({1, 0, 1}).normalized();
            domain.addBoundaryNode(point, RIGHT, normal);

            scalar_t inner_step = (end_[2] - beg_[2]) / (counts[2] - 1);  // z
            for (int j = 1; j < counts[2]-1; ++j) {
                point[2] = beg_[2] + j*inner_step;
                point[0] = beg_[0];
                normal = vec_t({-1, 0, 0});
                domain.addBoundaryNode(point, LEFT, normal);
                point[0] = end_[0];
                normal = vec_t({1, 0, 0});
                domain.addBoundaryNode(point, RIGHT, normal);
            }
        }
        step = (end_[2] - beg_[2]) / (counts[2] - 1);  // z
        for (int i = 1; i < counts[2]-1; ++i) {
            point[2] = beg_[2] + i*step;
            point[0] = beg_[0];
            point[1] = beg_[1];
            normal = vec_t({-1, -1, 0}).normalized();
            domain.addBoundaryNode(point, FRONT, normal);
            point[0] = beg_[0];
            point[1] = end_[1];
            normal = vec_t({-1, 1, 0}).normalized();
            domain.addBoundaryNode(point, LEFT, normal);
            point[0] = end_[0];
            point[1] = beg_[1];
            normal = vec_t({1, -1, 0}).normalized();
            domain.addBoundaryNode(point, RIGHT, normal);
            point[0] = end_[0];
            point[1] = end_[1];
            normal = vec_t({1, 1, 0}).normalized();
            domain.addBoundaryNode(point, BACK, normal);
        }
        domain.addBoundaryNode(beg_, BOTTOM, vec_t({-1, -1, -1}).normalized());
        domain.addBoundaryNode(end_, TOP, vec_t({1, 1, 1}).normalized());
        domain.addBoundaryNode(vec_t({end_[0], beg_[1], beg_[2]}), BOTTOM,
                               vec_t({1, -1, -1}).normalized());
        domain.addBoundaryNode(vec_t({beg_[0], end_[1], end_[2]}), TOP,
                               vec_t({-1, 1, 1}).normalized());
        domain.addBoundaryNode(vec_t({end_[0], end_[1], beg_[2]}), BACK,
                               vec_t({1, 1, -1}).normalized());
        domain.addBoundaryNode(vec_t({beg_[0], end_[1], beg_[2]}), LEFT,
                               vec_t({-1, 1, -1}).normalized());
        domain.addBoundaryNode(vec_t({beg_[0], beg_[1], end_[2]}), FRONT,
                               vec_t({-1, -1, 1}).normalized());
        domain.addBoundaryNode(vec_t({end_[0], beg_[1], end_[2]}), RIGHT,
                               vec_t({1, -1, 1}).normalized());
    }
    return domain;
}

template <typename vec_t>
std::ostream& BoxShape<vec_t>::print(std::ostream& os) const {
    return os << "BoxShape(" << beg_.transpose() << ", " << end_.transpose() << ")";
}

template <typename vec_t>
DomainDiscretization<vec_t>
BoxShape<vec_t>::discretizeBoundaryWithDensity(const std::function<scalar_t(vec_t)>& dr,
                                               int type) const {
    static_assert(0 <= dim && dim <= 3, "Only dimensions up to 3 are supported.");
    DomainDiscretization<vec_t> domain(*this);
    if (dim == 1) {
        domain.addBoundaryNode(beg_, type == 0 ? -1 : type, vec_t(-1.0));
        domain.addBoundaryNode(end_, type == 0 ? -2 : type, vec_t(1.0));
    } else if (dim == 2) {
        const int TOP = (type == 0) ? -4 : type;
        const int RIGHT = (type == 0) ? -2 : type;
        const int BOTTOM = (type == 0) ? -3 : type;
        const int LEFT = (type == 0) ? -1 : type;

        vec_t p1 = beg_;
        vec_t p2 = {beg_[0], end_[1]};
        vec_t p3 = end_;
        vec_t p4 = {end_[0], beg_[1]};

        for (const auto& p : discretization_helpers::discretizeLineWithDensity(p1, p2, dr)) {
            domain.addBoundaryNode(p, LEFT, {-1, 0});
        }
        for (const auto& p : discretization_helpers::discretizeLineWithDensity(p2, p3, dr)) {
            domain.addBoundaryNode(p, TOP, {0, 1});
        }
        for (const auto& p : discretization_helpers::discretizeLineWithDensity(p3, p4, dr)) {
            domain.addBoundaryNode(p, RIGHT, {1, 0});
        }
        for (const auto& p : discretization_helpers::discretizeLineWithDensity(p4, p1, dr)) {
            domain.addBoundaryNode(p, BOTTOM, {0, -1});
        }

        domain.addBoundaryNode(p1, LEFT, vec_t({-1, -1}).normalized());
        domain.addBoundaryNode(p2, TOP, vec_t({-1, 1}).normalized());
        domain.addBoundaryNode(p3, RIGHT, vec_t({1, 1}).normalized());
        domain.addBoundaryNode(p4, BOTTOM, vec_t({1, -1}).normalized());
    } else if (dim == 3) {
        int LEFT = -1, RIGHT = -2, FRONT = -3, BACK = -4, BOTTOM = -5, TOP = -6;
        if (type != 0) LEFT = RIGHT = FRONT = BACK = BOTTOM = TOP = type;
        Vec<scalar_t, 2> side_beg, side_end;
        GeneralFill<Vec<scalar_t, 2>> fill;
        fill.seed(0).proximityTolerance(0.99);
        {
            side_beg[0] = beg_[1]; side_beg[1] = beg_[2];
            side_end[0] = end_[1]; side_end[1] = end_[2];
            BoxShape<Vec<scalar_t, 2>> side(side_beg, side_end);
            auto dr2 = [&](const Vec<scalar_t, 2>& p) { return dr({beg_[0], p[0], p[1]}); };
            auto domain2d = side.discretizeWithDensity(dr2, fill);
            vec_t pos; pos[0] = beg_[0];
            for (int i : domain2d.interior()) {
                pos[1] = domain2d.pos(i, 0); pos[2] = domain2d.pos(i, 1);
                domain.addBoundaryNode(pos, LEFT, {-1, 0, 0});
            }
            auto dr22 = [&](const Vec<scalar_t, 2>& p) { return dr({end_[0], p[0], p[1]}); };
            domain2d = side.discretizeWithDensity(dr22, fill);
            pos[0] = end_[0];
            for (int i : domain2d.interior()) {
                pos[1] = domain2d.pos(i, 0); pos[2] = domain2d.pos(i, 1);
                domain.addBoundaryNode(pos, RIGHT, {1, 0, 0});
            }
        }
        {
            side_beg[0] = beg_[0]; side_beg[1] = beg_[2];
            side_end[0] = end_[0]; side_end[1] = end_[2];
            BoxShape<Vec<scalar_t, 2>> side(side_beg, side_end);
            auto dr2 = [&](const Vec<scalar_t, 2>& p) { return dr({p[0], beg_[1], p[1]}); };
            auto domain2d = side.discretizeWithDensity(dr2, fill);
            vec_t pos; pos[1] = beg_[1];
            for (int i : domain2d.interior()) {
                pos[0] = domain2d.pos(i, 0); pos[2] = domain2d.pos(i, 1);
                domain.addBoundaryNode(pos, FRONT, {0, -1, 0});
            }
            auto dr22 = [&](const Vec<scalar_t, 2>& p) { return dr({p[0], end_[1], p[1]}); };
            domain2d = side.discretizeWithDensity(dr22, fill);
            pos[1] = end_[1];
            for (int i : domain2d.interior()) {
                pos[0] = domain2d.pos(i, 0); pos[2] = domain2d.pos(i, 1);
                domain.addBoundaryNode(pos, BACK, {0, 1, 0});
            }
        }
        {
            side_beg[0] = beg_[0]; side_beg[1] = beg_[1];
            side_end[0] = end_[0]; side_end[1] = end_[1];
            BoxShape<Vec<scalar_t, 2>> side(side_beg, side_end);
            auto dr2 = [&](const Vec<scalar_t, 2>& p) { return dr({p[0], p[1], beg_[2]}); };
            auto domain2d = side.discretizeWithDensity(dr2, fill);
            vec_t pos; pos[2] = beg_[2];
            for (int i : domain2d.interior()) {
                pos[0] = domain2d.pos(i, 0); pos[1] = domain2d.pos(i, 1);
                domain.addBoundaryNode(pos, BOTTOM, {0, 0, -1});
            }
            auto dr22 = [&](const Vec<scalar_t, 2>& p) { return dr({p[0], p[1], end_[2]}); };
            domain2d = side.discretizeWithDensity(dr22, fill);
            pos[2] = end_[2];
            for (int i : domain2d.interior()) {
                pos[0] = domain2d.pos(i, 0); pos[1] = domain2d.pos(i, 1);
                domain.addBoundaryNode(pos, TOP, {0, 0, 1});
            }
        }
    }
    return domain;
}

template <typename vec_t>
DomainDiscretization<vec_t> BoxShape<vec_t>::discretizeWithStep(
        scalar_t step, int internal_type, int boundary_type) const {
    Vec<int, dim> counts;
    for (int i = 0; i < dim; ++i) counts[i] = iceil((end_[i] - beg_[i]) / step) + 1;
    return discretize(counts, internal_type, boundary_type);
}

template <typename vec_t>
DomainDiscretization<vec_t> BoxShape<vec_t>::discretize(
        const Vec<int, dim>& counts, int internal_type, int boundary_type) const {
    auto domain = discretizeBoundary(counts, boundary_type);
    Vec<int, dim> counts_internal = counts - Vec<int, dim>::Constant(2);
    if (internal_type == 0) internal_type = 1;
    for (const vec_t& p : linspace(beg_, end_, counts_internal, false)) {
        domain.addInternalNode(p, internal_type);
    }
    return domain;
}

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_BOXSHAPE_HPP_
