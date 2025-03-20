#ifndef MEDUSA_BITS_DOMAINS_BALLSHAPE_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_BALLSHAPE_FWD_HPP_

/**
 * @file
 * Declaration of class for ball shaped domains.
 *
 * @example test/domains/BallShape_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/types/Vec_fwd.hpp>
#include "DomainShape_fwd.hpp"

namespace mm {
/**
 * Class for working with ball shaped domains.
 *
 * Usage example:
 * @snippet domains/BallShape_test.cpp BallShape usage example
 * @ingroup domains
 */
template <typename vec_t>
class BallShape : public DomainShape<vec_t> {
    vec_t center_;  ///< Center of the ball.
    double radius_;  ///< Radius of the ball.
    using DomainShape<vec_t>::margin_;

  public:
    using typename DomainShape<vec_t>::scalar_t;
    using DomainShape<vec_t>::dim;
    using DomainShape<vec_t>::discretizeBoundaryWithDensity;
    using DomainShape<vec_t>::discretizeWithDensity;
    using DomainShape<vec_t>::discretizeBoundaryWithStep;
    using DomainShape<vec_t>::discretizeWithStep;
    /**
     * Constructs a `d`-dimensional ball defined by its center and radius.
     * @param center Position of the centre of the ball.
     * @param radius Radius of the ball.
     */
    BallShape(const vec_t& center, double radius) : center_(center), radius_(radius) {}
    /// Returns the position of the centre of the ball.
    const vec_t& center() const { return center_; }
    /// Returns the radius of the ball.
    scalar_t radius() const { return radius_; }

    bool contains(const vec_t& point) const override {
        return (center_ - point).squaredNorm() < (radius_ + margin_)*(radius_ + margin_);
    }

    std::pair<vec_t, vec_t> bbox() const override {
        return {center_ - vec_t(radius_), center_ + vec_t(radius_)};
    }

    DomainDiscretization<vec_t> discretizeBoundaryWithStep(scalar_t step, int type) const override;
    DomainDiscretization<vec_t> discretizeWithStep(
            scalar_t step, int internal_type, int boundary_type) const override;
    DomainDiscretization<vec_t> discretizeBoundaryWithDensity(
            const std::function<scalar_t(vec_t)>& dr, int type) const override;

    BallShape<vec_t>* clone() const override { return new BallShape<vec_t>(*this); }
    std::ostream& print(std::ostream& os) const override;
};

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_BALLSHAPE_FWD_HPP_
