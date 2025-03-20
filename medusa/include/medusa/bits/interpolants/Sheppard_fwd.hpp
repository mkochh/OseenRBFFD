#ifndef MEDUSA_BITS_INTERPOLANTS_SHEPPARD_FWD_HPP_
#define MEDUSA_BITS_INTERPOLANTS_SHEPPARD_FWD_HPP_

/**
 * @file
 * Declaration of class for Sheppard's interpolant.
 *
 * @example test/interpolants/SheppardInterpolant_test.cpp
 */

#include <medusa/bits/spatial_search/KDTree.hpp>
#include <medusa/bits/types/Range_fwd.hpp>

namespace mm {

/**
 * Scattered interpolant using a slightly modified Sheppard's interpolation (inverse distance
 * weighting).
 *
 * The essential difference is the introduction of a "regularization" parameter (`reg`) that is used
 * when computing the inverse distance weights:
 * \f$ w_i = 1 / (d(x, x_i)^{power} + reg) \f$
 * By supplying a positive regularization parameter one can avoid singularities at the locations of
 * the data points as well as control the "smoothness" of the interpolation (e.g., make the weights
 * of the neighbors less varied). The "smoothness" of interpolation can also be controlled by the
 * power parameter (`power`). The interpolation is not continuous due to the cut off at the closest
 * neighbors.
 *
 * @sa PUApproximant
 *
 * Usage example:
 * @snippet interpolants/SheppardInterpolant_test.cpp Sheppard interpolant usage example
 * @ingroup interpolants
 */
template <class vec_t, class value_t>
class SheppardInterpolant {
  private:
    KDTree<vec_t> tree;                         ///< Tree of all points.
    Range<value_t> values;                      ///< Function values at given points.
    typedef typename vec_t::scalar_t scalar_t;  ///< Scalar data type.

  public:
    /**
     * Construct a new empty Sheppard Interpolant object. Use @ref setValues and
     * @ref setPositions() to fill the object.
     */
    SheppardInterpolant() = default;

    /**
     * Construct a new Sheppard Interpolant. Use @ref setValues to add values.
     * @param pos Positions.
     */
    explicit SheppardInterpolant(const Range<vec_t>& pos) : tree(pos) {
        int N = pos.size();
        assert_msg(N > 0, "Number of positions in the tree must be greater than 0, got %d.", N);
    }

    /**
     * Construct a new Sheppard Interpolant object.
     * @param pos Initial positions.
     * @param values Function values to be interpolated..
     */
    SheppardInterpolant(const Range<vec_t>& pos, const Range<value_t>& values)
            : tree(pos), values(values) {
        int N = pos.size();
        int M = values.size();

        assert_msg(N > 0, "Number of positions in the tree must be greater than 0, got %d.", N);
        assert_msg(N == M,
                   "Number of positions must equal number of values. Got %d positions in the tree "
                   "and %d values.",
                   N, M);
    }

    /**
     * Set new values. Any container that supports forward iteration with `.begin()` and `.end()`
     * can be used. The contained element type must be compatible with `value_t`.
     */
    template <typename values_container_t>
    void setValues(const values_container_t& new_values) {
        values = Range<value_t>(new_values.begin(), new_values.end());
    }

    /// Set new positions.
    void setPositions(const Range<vec_t>& pos) { tree.reset(pos); }

    /**
     * Evaluate the interpolant at the given point.
     *
     * The complexity of the evaluation function is
     * \f$ \mathcal{O}(m\log{n}) \f$
     * for \f$ m \f$ closest nodes and \f$ n \f$ positions in the k-d tree.
     * 
     * The approximation order equals \f$ O(h) \f$.
     *
     * @param point Location where evaluation is needed.
     * @param num_closest Number of closest neighbors.
     * @param power The power of the inverse distance used for the interpolation weights.
     * @param reg Regularization parameter.
     * @param conf_dist The confusion distance below which the interpolator should use the value of
     * the closest data point instead of attempting to interpolate.
     * @return value_t Interpolation value.
     */
    value_t operator()(const vec_t& point, int num_closest, int power = 2, double reg = 0.0,
                       double conf_dist = 1e-12) const;
};

}  // namespace mm

#endif  // MEDUSA_BITS_INTERPOLANTS_SHEPPARD_FWD_HPP_
