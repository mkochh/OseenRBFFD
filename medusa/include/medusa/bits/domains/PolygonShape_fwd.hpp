#ifndef MEDUSA_BITS_DOMAINS_POLYGONSHAPE_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_POLYGONSHAPE_FWD_HPP_

#include <medusa/Config.hpp>
#include "DomainShape_fwd.hpp"

/**
 * @file
 * Declaration of PolygonShape class.
 *
 * @example test/domains/PolygonShape_test.cpp
 */

namespace mm {

/**
 * Shape representing a simple (i.e.\ non self intersecting) nonempty
 * polygon in 2D, which is given as a sequence of points.
 *
 * @tparam vec_t This shape class is used only for 2D domains.
 *
 * Usage example:
 * @snippet domains/PolygonShape_test.cpp PolygonShape usage example
 * @ingroup domains
 */
template <typename vec_t>
class PolygonShape : public DomainShape<vec_t> {
    static_assert(vec_t::dim == 2, "Only available in 2 dimensions.");
    using base_t = DomainShape<vec_t>;  ///< Base class type.

    std::vector<vec_t> points_;  ///< The points that define the polygon, stored in CCW order.
    /// Polygon extended by margin_, for contains checks. It updates along with margin and
    /// has the first point repeated at the end.
    std::vector<vec_t> points_with_margin_;
    using DomainShape<vec_t>::margin_;

  public:
    using typename base_t::scalar_t;
    using typename base_t::vector_t;
    using base_t::dim;
    using base_t::discretizeBoundaryWithDensity;
    using base_t::discretizeBoundaryWithStep;

    /**
     * Create polygon given a sequence of points. The points are stored in CCW order.
     * It is the user's responsibility to ensure the polygon is non self intersecting.
     * Polygon area is computed to ensure that polygon is nonempty and to possibly reverse the
     * point ordering, so that points are stored in counter clockwise (CCW) order.
     */
    explicit PolygonShape(const std::vector<vec_t>& points);

    /// Get points representing the polygon.
    const std::vector<vec_t>& points() const { return points_; }

    /// Computes the polygon extended by margin as well.
    void setMargin(scalar_t margin) override;

    /**
     * Winding number test for point in a polygon inclusion (paying respect to margin).
     * Loosely based on: http://geomalgorithms.com/a03-_inclusion.html
     */
    bool contains(const vec_t& point) const override;

    std::pair<vec_t, vec_t> bbox() const override;

    /**
     * @copydoc base_t::discretizeBoundaryWithStep
     *
     * The boundary points are guaranteed to be added in counter clock-wise order.
     */
    DomainDiscretization<vec_t> discretizeBoundaryWithStep(scalar_t step, int type) const override;

    /**
     * @copydoc base_t::discretizeBoundaryWithStep
     *
     * The boundary points are guaranteed to be added in counter clock-wise order.
     */
     DomainDiscretization<vec_t> discretizeBoundaryWithDensity(
            const std::function<scalar_t(vec_t)>& dr, int type) const override;

    PolygonShape<vec_t>* clone() const override { return new PolygonShape<vec_t>(*this); }

    std::ostream& print(std::ostream& os) const override {
        return os << "Polygon shape with " << points_.size() << " points.";
    }

  private:
     /**
      * Tests if a point is left, on, or right of an infinite line.
      * Input:  three points P0, P1, and P2
      * @returns a number > 0 if P2 is left of the line through P0 and P1 (looking from P0 to P1),
      * 0 for P2 on the line P0P1 or a number < 0 for P2 right of the line P0P1.
      */
    static inline scalar_t isLeft(const vec_t& P0, const vec_t& P1, const vec_t& P2) {
        return (P1[0] - P0[0]) * (P2[1] - P0[1]) - (P2[0] - P0[0]) * (P1[1] - P0[1]);
    }

    /// Save a version of points extended by margin_ for later contains checks.
    void extendByMargin();
};

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_POLYGONSHAPE_FWD_HPP_
