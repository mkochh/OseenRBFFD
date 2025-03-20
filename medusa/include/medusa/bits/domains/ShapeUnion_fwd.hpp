#ifndef MEDUSA_BITS_DOMAINS_SHAPEUNION_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_SHAPEUNION_FWD_HPP_

#include <medusa/Config.hpp>
#include "DomainShape_fwd.hpp"
#include <medusa/bits/types/Vec_fwd.hpp>

/**
 * @file
 * Declaration of ShapeUnion class.
 *
 * @example test/domains/ShapeUnion_test.cpp
 */

namespace mm {

/**
 * Class representing a union of two shapes. Used when constructing shapes
 * for later discretization or obtained implicitly with
 * DomainShape::operator+(), DomainShape::add() and DomainDiscretization::add().
 *
 * Usage example:
 * @snippet domains/ShapeUnion_test.cpp ShapeUnion usage example
 * @sa ShapeUnion, DomainShape, DomainDiscretization
 * @ingroup domains
 */
template <typename vec_t>
class ShapeUnion : public DomainShape<vec_t> {
    deep_copy_unique_ptr<DomainShape<vec_t>>
            sh1,  ///< First shape in the union.
            sh2;  ///< Second shape in the union.

  public:
    using typename DomainShape<vec_t>::scalar_t;
    using DomainShape<vec_t>::dim;
    using DomainShape<vec_t>::discretizeBoundaryWithDensity;
    using DomainShape<vec_t>::discretizeWithDensity;
    using DomainShape<vec_t>::discretizeBoundaryWithStep;
    using DomainShape<vec_t>::discretizeWithStep;

    /// Construct a union of two shapes.
    ShapeUnion(const DomainShape<vec_t>& sh1, const DomainShape<vec_t>& sh2) : sh1(sh1), sh2(sh2) {}

    /// Access the first shape.
    const DomainShape<vec_t>& first() const { return *sh1; }
    /// Access the second shape.
    const DomainShape<vec_t>& second() const { return *sh2; }

    bool contains(const vec_t& point) const override {
        return sh1->contains(point) || sh2->contains(point);
    }

    bool hasContains() const override {
        return sh1->hasContains() && sh2->hasContains();
    }

    std::pair<vec_t, vec_t> bbox() const override;

    DomainDiscretization<vec_t> discretizeBoundaryWithStep(scalar_t step, int type) const override;
    DomainDiscretization<vec_t> discretizeWithStep(
            scalar_t step, int internal_type, int boundary_type) const override;
    DomainDiscretization<vec_t> discretizeBoundaryWithDensity(
            const std::function<scalar_t(vec_t)>& dr, int type) const override;
    DomainDiscretization<vec_t> discretizeWithDensity(
            const std::function<scalar_t(vec_t)>& dr, int internal_type,
            int boundary_type) const override;

    ShapeUnion<vec_t>* clone() const override { return new ShapeUnion<vec_t>(*this); }
    std::ostream& print(std::ostream& os) const override;
};

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_SHAPEUNION_FWD_HPP_
