#ifndef MEDUSA_BITS_DOMAINS_SHAPEUNION_HPP_
#define MEDUSA_BITS_DOMAINS_SHAPEUNION_HPP_

#include "ShapeUnion_fwd.hpp"
#include "DomainDiscretization.hpp"
#include <utility>
#include <cmath>

/**
 * @file
 * Implementation of ShapeUnion class.
 */

namespace mm {

template <typename vec_t>
DomainDiscretization<vec_t> ShapeUnion<vec_t>::discretizeBoundaryWithStep(scalar_t step,
                                                                          int type) const {
    auto d1 = sh1->discretizeBoundaryWithStep(step, type);
    auto d2 = sh2->discretizeBoundaryWithStep(step, type);
    return d1.add(d2);
}

template <typename vec_t>
DomainDiscretization<vec_t> ShapeUnion<vec_t>::discretizeWithStep(
        scalar_t step, int internal_type, int boundary_type) const {
    auto d1 = sh1->discretizeWithStep(step, internal_type, boundary_type);
    auto d2 = sh2->discretizeWithStep(step, internal_type, boundary_type);
    return d1.add(d2);
}

template <typename vec_t>
DomainDiscretization<vec_t> ShapeUnion<vec_t>::discretizeBoundaryWithDensity(
        const std::function<scalar_t(vec_t)>& dr, int type) const {
    auto d1 = sh1->discretizeBoundaryWithDensity(dr, type);
    auto d2 = sh2->discretizeBoundaryWithDensity(dr, type);
    return d1.add(d2);
}

template <typename vec_t>
DomainDiscretization<vec_t> ShapeUnion<vec_t>::discretizeWithDensity(
        const std::function<scalar_t(vec_t)>& dr, int internal_type,
        int boundary_type) const {
    auto d1 = sh1->discretizeWithDensity(dr, internal_type, boundary_type);
    auto d2 = sh2->discretizeWithDensity(dr, internal_type, boundary_type);
    return d1.add(d2);
}

template <typename vec_t>
std::pair<vec_t, vec_t> ShapeUnion<vec_t>::bbox() const {
    vec_t lo, hi, lo1, hi1, lo2, hi2;
    std::tie(lo1, hi1) = sh1->bbox();
    std::tie(lo2, hi2) = sh2->bbox();
    for (int i = 0; i < dim; ++i) {
        lo[i] = std::min(lo1[i], lo2[i]);
        hi[i] = std::max(hi1[i], hi2[i]);
    }
    return {lo, hi};
}

template <typename vec_t>
std::ostream& ShapeUnion<vec_t>::print(std::ostream& os) const {
    return os << "ShapeUnion(" << *sh1 << ", " << *sh2 << ")";
}

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_SHAPEUNION_HPP_
