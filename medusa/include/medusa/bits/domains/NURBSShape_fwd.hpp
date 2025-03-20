#ifndef MEDUSA_BITS_DOMAINS_NURBSSHAPE_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_NURBSSHAPE_FWD_HPP_

/**
 * @file
 * Declaration of NURBSShape class.
 *
 * @example test/domains/NURBSShape_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/domains/DomainShape_fwd.hpp>
#include <medusa/bits/types/Vec.hpp>
#include <medusa/bits/types/Range.hpp>
#include <medusa/bits/domains/NURBSPatch_fwd.hpp>
#include <medusa/bits/utils/randutils.hpp>

namespace mm {

/**
 * Implementation details of NURBSShape, contains structures for partial class specialization.
 */
namespace nurbs_shape_internal {
/**
 * Internal structure of NURBSShape that helps with partial class specialization.
 * All of its functions have the same signature (apart from accepting a reference to NURBSShape)
 * as the corresponding functions of NURBSShape. NURBSShape holds an instance of this class
 * and calls its functions when specialized functions are requested.
 *
 * Note: partial class specialization can also be done by having a base class with functions,
 * that are not specialized (NURBSShapeBase) and a fully specialized class (NURBSShape) that
 * inherits from it. This approach was chosen to keep Medusa's class hierarchy cleaner.
 */
template <typename vec_t, typename param_vec_t>
struct NURBSShapeHelper {};
}  // namespace nurbs_shape_internal

/**
 * Class representing a shape made out of NURBS patches in an arbitrary dimensional space.
 *
 * @warning Currently only supports NURBS surfaces (2D domain) and NURBS
 * curves (1D domain).
 *
 * Usage example:
 * @snippet domains/NURBSShape_test.cpp NURBSShape usage example
 *
 * @tparam vec_t Vector type.
 *
 * @ingroup domains
 */
template <typename vec_t, typename param_vec_t>
class NURBSShape : public DomainShape<vec_t> {
    friend struct nurbs_shape_internal::NURBSShapeHelper<vec_t, param_vec_t>;

  public:
    typedef typename vec_t::scalar_t scalar_t;  ///< Scalar type.
    /// Store dimensionality.
    enum { /** Dimensionality of space after projection. */ dim = vec_t::dim,
        /** Dimensionality of space before projection. */ proj_dim = dim + 1,
        /** Dimensionality of parametric space. */ param_dim = param_vec_t::dim};
    typedef Vec<scalar_t, proj_dim> proj_vec_t;  ///< Vector type before projection.

  private:
    Range<NURBSPatch<vec_t, param_vec_t>> patches;  ///< Range of NURBS patches.
    int max_points = 5000000;  ///< Maximal number of points generated in surface fill.
    int seed_;  ///<  Seed for the random number generator.
    int n_samples = 15;  ///< Number of samples in surface fill.
    scalar_t zeta = 1 - 1e-10;  ///< Proximity tolerance in surface fill.
    scalar_t epsilon = 0;  ///< Evaluate normals slightly away from the boundary.

  public:
    /// Construct NURBSShape from a Range of patches.
    NURBSShape(const Range<NURBSPatch<vec_t, param_vec_t>>& patches_in);

    /// Move constructor for constructing NURBSShape from a Range of patches
    NURBSShape(Range<NURBSPatch<vec_t, param_vec_t>>&& patches_in);

    bool contains(const vec_t& /* point */) const override { return true; }

    bool hasContains() const override { return false; }

    std::pair<vec_t, vec_t> bbox() const override {
        assert_msg(false, "This function is not available for this shape. A"
                          "bounding box of its DomainDiscretization should be taken.");
        return {};
    }

    DomainDiscretization<vec_t> discretizeBoundaryWithDensity(
            const std::function<scalar_t(vec_t)>&, int) const override;

    std::ostream& print(std::ostream& ostream) const override {
        return ostream << "NURBS shape made of " << patches.size() << " patches.";
    }

    NURBSShape* clone() const override {
        return new NURBSShape<vec_t, param_vec_t>(*this);
    }

    /** Maximal number of points generated in surface fill
     * @warning The actual upper bound for number of points generated is maxPoints @f$\pm@f$
     * numSamples.
     */
    NURBSShape& maxPoints(int max_points) { this->max_points = max_points; return *this; }

    /// Set custom seed for the random number generator.
    NURBSShape& seed(int seed) { seed_ = seed; return *this; }

    /** Set proximity tolerance in surface fill.
     * A new candidate mush be at least `zeta*h(p)` away. It should hold that `0 < zeta < 1`.
     */
    NURBSShape& proximityTolerance(scalar_t zeta);

    /**
     * Controls the number of generated candidates from each point in surface fill.
     * For 2-D it is the actual number
     * of candidates, for 3-D it is the number of candidates on the great circle. Its value
     * is ignored in 1-D.
     */
    NURBSShape& numSamples(int n_samples) { this->n_samples = n_samples; return *this; }

    /**
     * Calculate boundary normals `epsilon` times the size of domain away from the boundary.
     * By default `epsilon = 0`. This can be useful if normals are undefined on the boundary.
     * It should hold that `0 < epsilon < 1`.
     */
    NURBSShape& boundaryProximity(scalar_t epsilon);
};


}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_NURBSSHAPE_FWD_HPP_
