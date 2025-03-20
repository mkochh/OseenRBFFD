#ifndef MEDUSA_BITS_DOMAINS_TRANSLATEDSHAPE_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_TRANSLATEDSHAPE_FWD_HPP_

#include <medusa/Config.hpp>
#include <medusa/bits/types/Vec_fwd.hpp>
#include "DomainShape_fwd.hpp"

/**
 * @file
 * Declaration of class for translated domain shapes.
 *
 * @example test/domains/TranslatedShape_test.cpp
 */

namespace mm {

/**
 * Class for working with translated domain shapes. An existing domain shape can be translated
 * to another point in the space, along with correctly working `contains` and discretization
 * methods.
 *
 * Usage example:
 * @snippet domains/TranslatedShape_test.cpp TranslatedShape usage example
 * @ingroup domains
 */
template <typename vec_t>
class TranslatedShape : public DomainShape<vec_t> {
    deep_copy_unique_ptr<DomainShape<vec_t>> sh;  ///< Shape to be translated.
    vec_t a;  ///< Translation vector

  public:
    /**
     * Construct a translated shape by specifying a shape and a translation vector.
     * Can be easily constructed using DomainShape::translate or DomainDiscretization::translate.
     */
    TranslatedShape(const DomainShape<vec_t>& sh, const vec_t& a);

    using typename DomainShape<vec_t>::scalar_t;
    using DomainShape<vec_t>::dim;
    using DomainShape<vec_t>::discretizeBoundaryWithDensity;
    using DomainShape<vec_t>::discretizeWithDensity;
    using DomainShape<vec_t>::discretizeBoundaryWithStep;
    using DomainShape<vec_t>::discretizeWithStep;

    bool contains(const vec_t& point) const override {
        return sh->contains(point-a);
    }

    bool hasContains() const override { return sh->hasContains(); }

    std::pair<vec_t, vec_t> bbox() const override {
        auto bb = sh->bbox();
        return {bb.first+a, bb.second+a};
    }

    TranslatedShape<vec_t>* clone() const override { return new TranslatedShape(*this); }

    std::ostream& print(std::ostream& os) const override {
        return os << "Shape " << *sh << " translated by " << a;
    }

    DomainDiscretization<vec_t> discretizeBoundaryWithStep(scalar_t step, int type) const override;

    DomainDiscretization<vec_t> discretizeWithStep(scalar_t step, int internal_type,
                                                   int boundary_type) const override;

    DomainDiscretization<vec_t> discretizeWithDensity(
            const std::function<scalar_t(vec_t)>& dr, int internal_type,
            int boundary_type) const override;

    DomainDiscretization<vec_t> discretizeBoundaryWithDensity(
            const std::function<scalar_t(vec_t)>& dr, int type) const override;

    /// Returns the underlying shape.
    const DomainShape<vec_t>& shape() const { return *sh; }

    /// Returns the translation vector.
    vec_t translation() const { return a; }
};

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_TRANSLATEDSHAPE_FWD_HPP_
