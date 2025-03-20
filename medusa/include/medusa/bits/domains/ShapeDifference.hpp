#ifndef MEDUSA_BITS_DOMAINS_SHAPEDIFFERENCE_HPP_
#define MEDUSA_BITS_DOMAINS_SHAPEDIFFERENCE_HPP_

#include "ShapeDifference_fwd.hpp"
#include "DomainDiscretization.hpp"
#include <ostream>
#include <medusa/bits/types/Vec.hpp>

/**
 * @file
 * Implementation of ShapeDifference class.
 */

namespace mm {

template <typename vec_t>
DomainDiscretization<vec_t> ShapeDifference<vec_t>::discretizeBoundaryWithStep(scalar_t step,
                                                                               int type) const {
    auto d1 = sh1->discretizeBoundaryWithStep(step, type);
    sh2->toggleMargin();
    auto d2 = sh2->discretizeBoundaryWithStep(step, type);
    sh2->toggleMargin();
    return d1.subtract(d2);
}

template <typename vec_t>
DomainDiscretization<vec_t> ShapeDifference<vec_t>::discretizeWithStep(
        scalar_t step, int internal_type, int boundary_type) const {
    auto d1 = sh1->discretizeWithStep(step, internal_type, boundary_type);
    sh2->toggleMargin();
    auto d2 = sh2->discretizeBoundaryWithStep(step, boundary_type);
    sh2->toggleMargin();
    return d1.subtract(d2);
}

template <typename vec_t>
DomainDiscretization<vec_t> ShapeDifference<vec_t>::discretizeBoundaryWithDensity(
        const std::function<scalar_t(vec_t)>& dr, int type) const {
    auto d1 = sh1->discretizeBoundaryWithDensity(dr, type);
    sh2->toggleMargin();
    auto d2 = sh2->discretizeBoundaryWithDensity(dr, type);
    sh2->toggleMargin();
    return d1.subtract(d2);
}

template <typename vec_t>
DomainDiscretization<vec_t> ShapeDifference<vec_t>::discretizeWithDensity(
        const std::function<scalar_t(vec_t)>& dr, int internal_type, int boundary_type) const {
    auto d1 = sh1->discretizeWithDensity(dr, internal_type, boundary_type);
    sh2->toggleMargin();
    auto d2 = sh2->discretizeBoundaryWithDensity(dr, boundary_type);
    sh2->toggleMargin();
    return d1.subtract(d2);
}

template <typename vec_t>
std::ostream& ShapeDifference<vec_t>::print(std::ostream& os) const {
    return os << "ShapeDifference(" << *sh1 << ", " << *sh2 << ")";
}

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_SHAPEDIFFERENCE_HPP_
