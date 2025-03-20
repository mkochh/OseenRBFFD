#ifndef MEDUSA_BITS_DOMAINS_DOMAINSHAPE_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_DOMAINSHAPE_FWD_HPP_

/**
 * @file
 * Declaration of base class for domain shapes.
 *
 * @example test/domains/DomainShape_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/utils/memutils.hpp>
#include <medusa/bits/types/Vec_fwd.hpp>
#include <iosfwd>
#include <utility>
#include <ostream>
#include <functional>

namespace mm {

template <class vec_t>
class deep_copy_unique_ptr;

template <class vec_t>
class DomainDiscretization;

// These must be hidden during a doxygen run due to a bug in processing forward declarations.
// See https://github.com/doxygen/doxygen/issues/8177
#ifndef DOXYGEN_RUN
template <typename vec_t>
class ShapeUnion;

template <typename vec_t>
class ShapeDifference;

template <typename vec_t>
class TranslatedShape;

template <typename vec_t>
class RotatedShape;
#endif

/**
 * Base class for geometric shapes of domains. This class implements the interface and
 * utilities for working with and discretizing geometric shapes.
 *
 * Usage example:
 * @snippet domains/DomainShape_test.cpp Domain shape usage example
 * @ingroup domains
 */
template <typename vec_t>
class DomainShape {
  public:
    typedef vec_t vector_t;  ///< Vector data type used in computations.
    typedef typename vec_t::Scalar scalar_t;   ///< Scalar data type used in computation.
    /// Store dimension of the domain.
    enum { /** Dimensionality of the domain. */ dim = vec_t::dim };

  protected:
    /**
     * Tolerance for the geometric operation of the domain. The domain should behave as if
     * it was `margin_` thicker. Default margin is `1e-10`.
     */
    scalar_t margin_;
  public:
    /// Construct domain with default margin.
    DomainShape() : margin_(1e-10) {}
    /// Virtual destructor to properly destruct base class when invoked polymorphically.
    virtual ~DomainShape() = default;

    /// Returns current margin.
    scalar_t margin() const { return margin_; }

    /// Sets domain margin to `margin`.
    virtual void setMargin(scalar_t margin) { margin_ = margin; }
    /// Toggles the margin from positive to negative.
    void toggleMargin() { setMargin(-margin()); }
    /// Returns a shape representing a union of `*this` and `other`.
    ShapeUnion<vec_t> add(const DomainShape& other) const;
    /// Operator form of DomainShape::add. @sa add
    ShapeUnion<vec_t> operator+(const DomainShape& other) const { return add(other); }

    /// Returns a shape representing a difference of `*this` and `other`.
    ShapeDifference<vec_t> subtract(const DomainShape& other) const;
    /// Operator form of DomainShape::subtract. @sa subtract
    ShapeDifference<vec_t> operator-(const DomainShape& other) const { return subtract(other); }

    /// Project point to boundary using bisection along the line define by `unit_normal`.
    virtual std::pair<bool, vec_t> projectPointToBoundary(const vec_t& point,
                                                          const vec_t& unit_normal) const;

    /// Return true if `point` is not more than `margin()` outside the domain.
    virtual bool contains(const vec_t& point) const = 0;

    /// Return true if shape has `contains()` method implemented
    virtual bool hasContains() const { return true; }

    /**
     * Return the bounding box of the domain. Bounding box is returned in format
     * `bbox() == {{mx, my, ...}, {MX, MY, ...}}`, such that `mx <= Mx` and `my <= My` etc.\
     * and that the whole domain is contained in the cuboid `[mx, my, ...] x [Mx, My, ...]`.
     */
    virtual std::pair<vec_t, vec_t> bbox() const = 0;
    /// Polymorphic clone pattern.
    virtual DomainShape* clone() const = 0;
    /// Output information about this shape to given output stream `os`.
    virtual std::ostream& print(std::ostream& os) const = 0;

    // Two types of discretizations: uniform step and density (will possibly add random later).

    // STEP -- calls density by default
    /**
     * Returns a discretization of the boundary of this shape with approximately uniform
     * distance `step` between nodes. `step` must be positive. Added nodes
     * are of type `type`, which must be non-positive. Value 0 indicates that types are
     * dependant on the implementation of concrete shape.
     */
    virtual DomainDiscretization<vec_t> discretizeBoundaryWithStep(scalar_t step, int type) const {
        return discretizeBoundaryWithDensity([=](const vec_t&) { return step; }, type);
    }

    /**
     * Returns a discretization of the boundary of this shape with approximately uniform
     * distance `step` between nodes. Node type are decided by the underlying shape.
     */
    DomainDiscretization<vec_t> discretizeBoundaryWithStep(scalar_t step) const {
        return this->discretizeBoundaryWithStep(step, 0);
    }
    /**
     * Returns a discretization of this shape with approximately uniform distance `step` between
     * nodes. `step` must be positive.
     * @param step Desired internodal distance.
     * @param internal_type User supplied type of internal nodes. Must be non-negative.
     * @param boundary_type User supplied type of boundary nodes. Must be non-positive.
     * If any of the types is 0, the underlying shape provides the default.
     */
    virtual DomainDiscretization<vec_t> discretizeWithStep(
            scalar_t step, int internal_type, int boundary_type) const {
        return discretizeWithDensity([=](const vec_t&) { return step; }, internal_type,
                                     boundary_type);
    }

    /// @ref discretizeWithStep but with default types as assigned by the shape.
    DomainDiscretization<vec_t> discretizeWithStep(scalar_t step) const {
        return this->discretizeWithStep(step, 0, 0);
    }

    // DENSITY -- calls fill engine in interior
    /**
     * Returns a discretization of the domain with spatially variable step.
     * @param dr Function giving desired internodal distance at each point.
     * @param internal_type User supplied type of internal nodes. Must be non-negative.
     * @param boundary_type User supplied type of boundary nodes. Must be non-positive.
     * If any of the types is 0, the underlying shape provides the default.
     * @return Discretization with nodes distributed according to `dr`.
     */
    virtual DomainDiscretization<vec_t> discretizeWithDensity(
            const std::function<scalar_t(vec_t)>& dr, int internal_type, int boundary_type) const;

    /// Overload for fill engine.
    template <typename func_t, typename fill_t>
    DomainDiscretization<vec_t> discretizeWithDensity(const func_t& dr, const fill_t& fill,
                                                      int internal_type, int boundary_type) const {
        DomainDiscretization<vec_t> domain = discretizeBoundaryWithDensity(dr, boundary_type);
        fill(domain, dr, internal_type);
        return domain;
    }

    /// Overload with default types.
    DomainDiscretization<vec_t> discretizeWithDensity(
            const std::function<scalar_t(vec_t)>& dr) const {
        return discretizeWithDensity(dr, 0, 0);
    }

    /// Overload for fill engine with default types.
    template <typename func_t, typename fill_t>
    DomainDiscretization<vec_t> discretizeWithDensity(const func_t& dr, const fill_t& fill) const {
        return discretizeWithDensity(dr, fill, 0, 0);
    }

    /**
     * Discretizes boundary with given density and fill engine.
     * If type is 0, the underlying shape provides the default.
     */
    virtual DomainDiscretization<vec_t> discretizeBoundaryWithDensity(
            const std::function<scalar_t(vec_t)>& dr, int type) const = 0;

    /// Overload with default type.
    DomainDiscretization<vec_t> discretizeBoundaryWithDensity(
            const std::function<scalar_t(vec_t)>& dr) const {
        return this->discretizeBoundaryWithDensity(dr, 0);
    }

    /**
     * Translate the shape by given vector `a`.
     * @note It is usually faster to first discretize the domain and than translate the whole
     * discretization using DomainDiscretization::translate.
     */
    TranslatedShape<vec_t> translate(const vec_t& a);

    /**
     * Transform the shape by given orthogonal matrix `Q`.
     * @note It is usually faster to first discretize the domain and than transform the whole
     * discretization using DomainDiscretization::rotate.
     */
    RotatedShape<vec_t> rotate(const Eigen::Matrix<scalar_t, dim, dim>& Q);

    /// 2D version of @ref rotate accepting an angle.
    RotatedShape<vec_t> rotate(scalar_t angle);

    template <typename V>
    friend std::ostream& operator<<(std::ostream& os, const DomainShape<V>& shape);
};

/**
 * Output info about given shape to ostream. Calls polymorphic DomainShape::print function.
 * @sa print
 */
template <typename V>
std::ostream& operator<<(std::ostream& os, const DomainShape<V>& shape) {
    return shape.print(os);
}

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_DOMAINSHAPE_FWD_HPP_
