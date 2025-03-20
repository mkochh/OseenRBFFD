#ifndef MEDUSA_BITS_DOMAINS_BOXSHAPE_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_BOXSHAPE_FWD_HPP_

/**
 * @file
 * Declaration of class for box shaped domains.
 *
 * @example test/domains/BoxShape_test.cpp
 */

#include <medusa/Config.hpp>
#include "DomainShape_fwd.hpp"
#include <medusa/bits/types/Range_fwd.hpp>
#include <medusa/bits/types/Vec_fwd.hpp>

namespace mm {

/**
 * Class for working with box shaped domains. This subclass of DomainShape implements the interface
 * for working with geometric box shaped domains (generalisation of a rectangle in arbitrary
 * number of dimensions).
 *
 * Usage example:
 * @snippet domains/BoxShape_test.cpp BoxShape usage example
 * @ingroup domains
 */
template <typename vec_t>
class BoxShape : public DomainShape<vec_t> {
    vec_t beg_,  ///< Smaller of the edge points.
          end_;  ///< Larger of the edge points.
    using typename DomainShape<vec_t>::scalar_t;
    using DomainShape<vec_t>::margin_;

  public:
    using DomainShape<vec_t>::dim;
    using DomainShape<vec_t>::discretizeBoundaryWithDensity;
    using DomainShape<vec_t>::discretizeWithDensity;
    using DomainShape<vec_t>::discretizeBoundaryWithStep;
    using DomainShape<vec_t>::discretizeWithStep;
    /**
     * Constructs a n dimensional box shaped domain defined by two n dimensional points -
     * vertices at ends of the longest diagonal.
     * @param beg Position vector of first point defining the box. (Vertex at lower end of longest
     * diagonal with every coordinate lower than `end` vertex).
     * @param end Position vector of second point defining the box. (Vertex at higher end of longest
     * diagonal with every coordinate higher than `beg` vertex).
     */
    BoxShape(vec_t beg, vec_t end);
    /// Returns the first point defining the box.
    const vec_t& beg() const { return beg_; }
    /// Returns the second point defining the box.
    const vec_t& end() const { return end_; }

    bool contains(const vec_t& point) const override;

    std::pair<vec_t, vec_t> bbox() const override { return {beg_, end_}; }

    DomainDiscretization<vec_t> discretizeBoundaryWithStep(scalar_t step, int type) const override;
    DomainDiscretization<vec_t> discretizeWithStep(
            scalar_t step, int internal_type, int boundary_type) const override;
    DomainDiscretization<vec_t> discretizeBoundaryWithDensity(
            const std::function<scalar_t(vec_t)>& dr, int type) const override;

    /**
     * Uniformly discretizes underlying domain boundary. It fills positions
     * with appropriate coordinates and types with distinctive negative integers
     * for each boundary of the box. Unless type is user specified,
     * all points where first coordinate (x) is constant have types -1 and -2, points
     * where second coordinate (y) is constant have types -3 and -4 and so on.
     * How types are assigned by default in 2D and 3D is presented on the included graphic.
     * @image html cube_rectangle.png
     * @image latex cube_rectangle.png
     *
     * @param counts Vector of same dimension as the domain.
     * Specifies how many discretization points to use in each dimension.
     * @param type User supplied type of boundary nodes. Must be non-positive.
     * If any types is 0, the default types as shown above are used.
     */
    DomainDiscretization<vec_t> discretizeBoundary(const Vec<int, dim>& counts, int type = 0) const;
    /**
     * Uniformly discretizes underlying domain. It fills positions
     * with appropriate coordinates and given type.
     * @param counts Vector of same dimension as the domain. Specifies how many discretization
     * points to use in each dimension.
     * @param internal_type User supplied type of internal nodes. Must be non-negative.
     * @param boundary_type User supplied type of boundary nodes. Must be non-positive.
     * If any of the types is 0, the underlying shape provides the default.
     */
    DomainDiscretization<vec_t> discretize(const Vec<int, dim>& counts, int internal_type = 0,
            int boundary_type = 0) const;

    BoxShape<vec_t>* clone() const override { return new BoxShape<vec_t>(*this); }

    std::ostream& print(std::ostream& os) const override;
};

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_BOXSHAPE_FWD_HPP_
