#ifndef MEDUSA_BITS_DOMAINS_TRANSLATEDSHAPE_HPP_
#define MEDUSA_BITS_DOMAINS_TRANSLATEDSHAPE_HPP_

/**
 * @file
 * Implementation of translated domain shapes.
 */

#include "TranslatedShape_fwd.hpp"
#include "DomainShape.hpp"
#include "DomainDiscretization.hpp"
#include <medusa/bits/utils/assert.hpp>

namespace mm {

template <typename vec_t>
TranslatedShape<vec_t>::TranslatedShape(const DomainShape<vec_t>& sh, const vec_t& a) :
        sh(sh), a(a) {
    auto* tsh = dynamic_cast<const TranslatedShape<vec_t>*>(&sh);
    if (tsh != nullptr) {  // collapse double translations
        this->sh = tsh->sh;
        this->a += tsh->a;
    }
}

template <typename vec_t>
DomainDiscretization <vec_t>
TranslatedShape<vec_t>::discretizeBoundaryWithStep(scalar_t step, int type) const  {
    auto d = sh->discretizeBoundaryWithStep(step, type);
    d.translate(a);
    return d;
}

template <typename vec_t>
DomainDiscretization<vec_t>
TranslatedShape<vec_t>::discretizeWithStep(scalar_t step, int internal_type,
                                           int boundary_type) const {
    auto d = sh->discretizeWithStep(step, internal_type, boundary_type);
    d.translate(a);
    return d;
}

template <typename vec_t>
DomainDiscretization<vec_t>
TranslatedShape<vec_t>::discretizeWithDensity(const std::function<scalar_t(vec_t)>& dr,
                                              int internal_type, int boundary_type) const {
    auto d = sh->discretizeWithDensity(dr, internal_type, boundary_type);
    d.translate(a);
    return d;
}

template <typename vec_t>
DomainDiscretization<vec_t>
TranslatedShape<vec_t>::discretizeBoundaryWithDensity(const std::function<scalar_t(vec_t)>& dr,
                                                      int type) const {
    auto d = sh->discretizeBoundaryWithDensity(dr, type);
    d.translate(a);
    return d;
}

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_TRANSLATEDSHAPE_HPP_
