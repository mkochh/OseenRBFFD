#ifndef MEDUSA_BITS_DOMAINS_DOMAINSHAPE_HPP_
#define MEDUSA_BITS_DOMAINS_DOMAINSHAPE_HPP_

/**
 * @file
 * Implementations of base class for domain shapes.
 */

#include "DomainShape_fwd.hpp"
#include "ShapeDifference.hpp"
#include "ShapeUnion.hpp"
#include "TranslatedShape.hpp"
#include "RotatedShape.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/utils/numutils.hpp>
#include <cmath>

namespace mm {

template <typename vec_t>
std::pair<bool, vec_t> DomainShape<vec_t>::projectPointToBoundary(
        const vec_t& point, const vec_t& unit_normal) const {
    assert_msg(unit_normal.norm() > 1e-9, "Normal %s given for point %s is zero.",
               unit_normal, point);
    if (dim == 1) return {false, point};
    // Find point on the other side of the boundary with exponential search
    vec_t start = point;
    vec_t normal = unit_normal;
    bool is_inside = contains(start);
    scalar_t max_stretch = 100 * margin_;
    while (true) {
        if (contains(start + max_stretch * normal) != is_inside) break;
        if (contains(start - max_stretch * normal) != is_inside) {
            max_stretch *= -1;
            break;
        }
        max_stretch *= 2;
        if (std::isinf(max_stretch)) {  // hint is bad
            return {false, vec_t()};
        }
    }

    // Find the point on the boundary using bisection
    scalar_t stretch = max_stretch;
    while (std::abs(stretch) > margin_) {
        stretch /= 2;
        if (contains(start + stretch * normal) == is_inside) {
            start += stretch * normal;
        }
    }

    // Make unit normal point outside
    if (is_inside) {
        normal *= signum(max_stretch);
    } else {
        normal *= -signum(max_stretch);
    }

    // Make sure point is inside
    while (!contains(start)) start -= margin_ * normal;

    return {true, start};
}

template <typename vec_t>
ShapeUnion<vec_t> DomainShape<vec_t>::add(const DomainShape& other) const {
    return ShapeUnion<vec_t>(*this, other);
}

template <typename vec_t>
ShapeDifference<vec_t> DomainShape<vec_t>::subtract(const DomainShape& other) const {
    return ShapeDifference<vec_t>(*this, other);
}

template <typename vec_t>
DomainDiscretization<vec_t>
DomainShape<vec_t>::discretizeWithDensity(const std::function<scalar_t(vec_t)>&, int, int) const {
    assert_msg(false, "This domain does not support filling with density, "
                      "use a general fill engine instead.");
    return DomainDiscretization<vec_t>(*this);
}

template <typename vec_t>
TranslatedShape<vec_t> DomainShape<vec_t>::translate(const vec_t& a) {
    return TranslatedShape<vec_t>(*this, a);
}

template <typename vec_t>
RotatedShape<vec_t> DomainShape<vec_t>::rotate(scalar_t angle) {
    assert_msg(dim == 2, "Angle rotation only available in 2D.");
    scalar_t s = std::sin(angle);
    scalar_t c = std::cos(angle);
    Eigen::Matrix<double, dim, dim> Q; Q << c, -s, s, c;
    return rotate(Q);
}

template <typename vec_t>
RotatedShape<vec_t> DomainShape<vec_t>::rotate(const Eigen::Matrix<scalar_t, dim, dim>& Q) {
    return RotatedShape<vec_t>(*this, Q);
}


}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_DOMAINSHAPE_HPP_
