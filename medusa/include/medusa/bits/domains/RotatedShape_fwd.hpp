#ifndef MEDUSA_BITS_DOMAINS_ROTATEDSHAPE_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_ROTATEDSHAPE_FWD_HPP_

/**
 * @file
 * Declaration of class for rotated domain shapes.
 *
 * @example test/domains/RotatedShape_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/types/Vec_fwd.hpp>
#include "DomainShape_fwd.hpp"

namespace mm {

/**
 * Class for working with rotated (or mirrored) domain shapes. An existing domain shape can be
 * transformed by any linear isometry (represented by an orthogonal matrix Q), along with
 * correctly working `contains` and discretization methods.
 *
 * Usage example:
 * @snippet domains/RotatedShape_test.cpp RotatedShape usage example
 * @ingroup domains
 */
template <typename vec_t>
class RotatedShape : public DomainShape<vec_t> {
  public:
    using typename DomainShape<vec_t>::scalar_t;
    using DomainShape<vec_t>::dim;

  private:
    deep_copy_unique_ptr<DomainShape<vec_t>> sh;  ///< Shape to be transformed.
    Eigen::Matrix<scalar_t, dim, dim> Q;  ///< Orthogonal transformation matrix.

  public:
    /**
     * Construct a transformed shape by specifying a shape and an orthogonal transformation matrix.
     * Can be easily constructed using DomainShape::rotate or DomainDiscretization::rotate.
     */
    RotatedShape(const DomainShape<vec_t>& sh, const Eigen::Matrix<scalar_t, dim, dim>& Q);

    using DomainShape<vec_t>::discretizeBoundaryWithDensity;
    using DomainShape<vec_t>::discretizeWithDensity;
    using DomainShape<vec_t>::discretizeBoundaryWithStep;
    using DomainShape<vec_t>::discretizeWithStep;

    /// Returns the underlying shape.
    const DomainShape<vec_t>& shape() const { return *sh; }

    /// Returns the transformation used in this shape.
    Eigen::Matrix<scalar_t, dim, dim> rotation() const { return Q; }

    bool contains(const vec_t& point) const override {
        return sh->contains(Q.transpose()*point);
    }

    bool hasContains() const override { return sh->hasContains(); }

    std::pair<vec_t, vec_t> bbox() const override;

    RotatedShape<vec_t>* clone() const override { return new RotatedShape(*this); }

    std::ostream& print(std::ostream& os) const override {
        return os << "Shape " << *sh << " rotated by " << Q;
    }

    /// @cond Doxgen parses this wrong, don't know why...
    DomainDiscretization<vec_t> discretizeBoundaryWithStep(scalar_t step, int type) const override;

    DomainDiscretization<vec_t> discretizeWithStep(scalar_t step, int internal_type,
                                                   int boundary_type) const override;

    DomainDiscretization<vec_t> discretizeWithDensity(
            const std::function<scalar_t(vec_t)>& dr, int internal_type,
            int boundary_type) const override;

    DomainDiscretization<vec_t> discretizeBoundaryWithDensity(
            const std::function<scalar_t(vec_t)>& dr, int type) const override;
    /// @endcond
};

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_ROTATEDSHAPE_FWD_HPP_
