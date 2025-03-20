#ifndef MEDUSA_BITS_INTERPOLANTS_PUAPPROXIMANT_FWD_HPP_
#define MEDUSA_BITS_INTERPOLANTS_PUAPPROXIMANT_FWD_HPP_

/**
 * @file
 * Declaration of class for Partition of Unity Approximant.
 *
 * @example test/interpolants/PUApproximant_test.cpp
 */

#include <Eigen/Core>
#include <medusa/bits/domains/DomainDiscretization_fwd.hpp>
#include <medusa/bits/types/Range_fwd.hpp>

namespace mm {

/**
 * An efficient partition-of-unity based approximation
 * method for gluing the local approximations together into a smooth field.
 *
 * To join the many different local approximations together, a weighted average of
 * the relevant solutions at each point is computed. This can be done continuously
 * and efficiently using a compactly supported partition of unity.
 *
 * There are more reasonable choices for \f$ r \f$, such as the distance to
 * the closest computation node, not included in the neighborhood. The larger
 * the radii \f$ r \f$ the more prominent the averaging effect will be.
 *
 * To evaluate approximant at \f$ N_q \f$ points, we can construct a spatial
 * search structure that supports radius based nearest-neighbor queries
 * (such as a \f$ k \f$-d tree).
 *
 * For more details see <a href="https://arxiv.org/abs/1610.07050">this paper</a>.
 *
 * @sa SheppardInterpolant
 *
 * @ingroup interpolants
 */
template <typename vec_t>
class PUApproximant {
    template <typename scalar_t>
    using eigen_vt = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;  ///< Vector type.
    using d_scalar_t = typename vec_t::scalar_t;                  ///< Scalar type.

  public:
    /**
     * Computes \f$ C^2 \f$ smooth weights \f$ w \in \left [ 0, 1 \right ] \f$.
     *
     * There are more reasonable choices for \f$ r \f$, such as the distance to
     * the closest computation node, not included in the neighborhood. The larger
     * the radii \f$ r \f$ the more prominent the averaging effect will be.
     *
     * @param r Effective radius.
     * @return vec_t::scalar_t
     */
    static typename vec_t::scalar_t weight(typename vec_t::scalar_t r) {
        assert_msg(r >= 0 && r <= 1,
                   "Effective radius must be larger than 0 and smaller than 1, got %g", r);
        return ipow<3>(1 - r) * (6 * r * r + 3 * r + 1);
    }

    /**
     * Evaluates the Partion-of-unity approximant at points 'query_points'
     * constructed on given \f$ (domain, values) \f$ pairs.
     *
     * The order of the approximation is equal to the order of local approximations
     * of \f$ engine \f$.
     *
     * Computational complexity is \f$ O((N + N_q) \log N_q) \f$, where \f$ N \f$ is the number of
     * points in the domain and \f$ N_q \f$ is the number of query points.
     *
     * Usage example:
     * @snippet interpolants/PUApproximant_test.cpp Partition-of-unity approximant usage example
     *
     * @param domain Domain with nodes and stencils.
     * @param values Function values in each domain node.
     * @param query_points Evaluation points where the approximated function
     * should be evaluated.
     * @param radius_factor Local subdomains have radius `radius_factor*h`, where `h`
     * is the distance to the farthest support node. Radius factor must be large enough that local
     * subdomains cover all query points.
     * @param engine Local approximation engine. If the engine is interpolatory, so is the
     * PU approximant.
     * @return Approximated function values in query points.
     */
    template <typename scalar_t, typename engine_t>
    static eigen_vt<scalar_t> evaluate(const DomainDiscretization<vec_t>& domain,
                                       const eigen_vt<scalar_t>& values,
                                       const Range<vec_t>& query_points, d_scalar_t radius_factor,
                                       const engine_t& engine);

    /// Overload of @ref evaluate using a precomputed \f$ k \f$-d tree.
    template <typename scalar_t, typename engine_t>
    static eigen_vt<scalar_t> evaluate(const DomainDiscretization<vec_t>& domain,
                                       const eigen_vt<scalar_t>& values,
                                       const KDTree<vec_t>& query_points_tree,
                                       d_scalar_t radius_factor, const engine_t& engine);
};

}  // namespace mm

#endif  // MEDUSA_BITS_INTERPOLANTS_PUAPPROXIMANT_FWD_HPP_
