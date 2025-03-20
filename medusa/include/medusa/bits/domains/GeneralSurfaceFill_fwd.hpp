#ifndef MEDUSA_BITS_DOMAINS_GENERALSURFACEFILL_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_GENERALSURFACEFILL_FWD_HPP_

/**
 * @file
 * Declaration of the general node placing algorithm for parametrically given surfaces.
 *
 * @example test/domains/GeneralSurfaceFill_test.cpp
 */

#include <medusa/Config.hpp>
#include <functional>

namespace mm {

template <typename vec_t>
class DomainDiscretization;

template <typename vec_t>
class DomainShape;

/**
 * Implements general `n`-d node placing algorithm for parametrically given `d`-d surfaces,
 * as described in https://arxiv.org/abs/2005.08767.
 * If you specifically use this algorithm, we would appreciate if you cite the above publication.
 *
 * Given a regular parametrization surface
 * @f$ \vec r: \Lambda \subseteq \mathbb{R}^d \to \partial \Omega \subseteq \mathbb{R}^n@f$,
 * and a spacing function @f$ h: \partial \Omega \to (0, \infty)@f$
 * the algorithm fills the domain discretization with points spaced approximately @f$h@f$ apart.
 *
 * The algorithm starts from existing nodes in the parametrization function's domain,
 * usually one random node. They are processed in a queue and for each node, new candidates are
 * generated around it on a circle. The circle's radius is determined in way that when the new
 * candidates are mapped by the parametrization function, their distance from the initial node is
 * approximately equal to the supplied density function. This is done using the first degree Taylor
 * polynomial. Candidates are then checked and some are accepted and added to the queue. This is
 * repeated until no more nodes can be generated or the maximal number of points has been reached.
 *
 * Usage example (filling the surface of a torus):
 * @snippet domains/GeneralSurfaceFill_test.cpp GeneralSurfaceFill 3d function fill usage example
 *
 * @warning If GeneralFill is not working after filling the boundary, normals
 * might be turned on the inside; see @ref compute_normal for more information.
 *
 * @ingroup domains
 */
template <typename vec_t, typename param_vec_t>
class GeneralSurfaceFill{
  public:
    typedef DomainDiscretization<vec_t> domain_t;  ///< Domain discretization type.
    /// Parametric domain discretization type.
    typedef DomainDiscretization<param_vec_t> param_domain_t;
    typedef typename vec_t::scalar_t scalar_t;  ///< Scalar type.
    typedef typename param_vec_t::scalar_t param_scalar_t;  ///< Parametric domain scalar type.
    /// Store dimension of the domains.
    enum { /** Dimensionality of the domain. */ dim = vec_t::dim,
        /** Dimensionality of the parametric domain. */ param_dim = param_vec_t::dim };

  private:
    int max_points = 5000000;  ///< Maximal number of points generated.
    int seed_;  ///<  Seed for the random number generator.
    int n_samples = 15;  ///< Number of samples.
    scalar_t zeta = 1 - 1e-10;  ///< Proximity tolerance.

  public:
    GeneralSurfaceFill();

    /** Maximal number of points generated
     * @warning The actual upper bound for number of points generated is maxPoints @f$\pm@f$
     * numSamples.
     */
    GeneralSurfaceFill& maxPoints(int max_points) { this->max_points = max_points; return *this; }
    /// Set custom seed for the random number generator.
    GeneralSurfaceFill& seed(int seed) { seed_ = seed; return *this; }
    /// Set proximity tolerance. A new candidate mush be at least `zeta*h(p)` away. It should hold
    /// that `0 < zeta < 1`.
    GeneralSurfaceFill& proximityTolerance(scalar_t zeta);
    /**
     * Controls the number of generated candidates from each point. For 2-D it is the actual number
     * of candidates, for 3-D it is the number of candidates on the great circle. Its value
     * is ignored in 1-D.
     */
    GeneralSurfaceFill& numSamples(int n_samples) { this->n_samples = n_samples; return *this; }

    /**
     * Fills the parametrically given surface with a quality node distribution according to
     * the `spacing_function`.
     * @param domain Domain @f$\Omega@f$ to fill with nodes. Can be partially filled.
     * @param param_domain Domain @f$\Lambda@f$ of the parametrization function. Can be partially
     * filled with parameters, which are taken as seeds. It is assumed that all parameters
     * have a corresponding node in `domain` (but domain can contain other nodes as well).
     * @param param_function Parametrization function @f$\vec r@f$.
     * @param param_jacobian Jacobian matrix @f$\nabla \vec r@f$ of the parametrization function.
     * Must have full rank.
     * @param spacing_function Nodal spacing function @f$h@f$.
     * @param tree Search structure on domain. This can be used if the domain already has
     * some nodes, such as when filling different patches separately.
     * @param type Type of the nodes. If 0, the engines default value is used.
     *
     * @warning Header <code> \#include <medusa/bits/domains/GeneralSurfaceFill.hpp></code> must be
     * included additionally to <code> \#include <medusa/Medusa_fwd.hpp></code>.
     */
    template <typename param_func_t, typename jm_func_t, typename search_t, typename spacing_func_t>
    void operator()(domain_t& domain, param_domain_t& param_domain,
                    const param_func_t& param_function, const jm_func_t& param_jacobian,
                    const spacing_func_t& spacing_function,
                    search_t& tree, int type = 0) const;

    /// Version with single function returning a pair of point and jacobian.
    template <typename param_func_t, typename search_t, typename spacing_func_t>
    void fillParametrization(domain_t& domain, param_domain_t& param_domain,
                    const param_func_t& param, const spacing_func_t& spacing_function,
                    search_t& tree, int type = 0) const;

    /// Fills given surface according to the nodal spacing function `spacing_function`.
    template <typename param_func_t, typename jm_func_t, typename spacing_func_t>
    void operator()(domain_t& domain, param_domain_t& param_domain,
                    const param_func_t& param_function, const jm_func_t& param_jacobian,
                    const spacing_func_t& spacing_function, int type = 0) const;

    /// Overload for constant function
    template <typename param_func_t, typename jm_func_t>
    void operator()(domain_t& domain, param_domain_t& param_domain,
                    const param_func_t& param_function, const jm_func_t& param_jacobian,
                    const scalar_t& h, int type = 0) const {
        this->operator()(domain, param_domain, param_function, param_jacobian,
                [=] (const vec_t&) { return h; }, type);
    }

    /// Overload for Shape instead of DomainDiscretization
    template <typename param_func_t, typename jm_func_t, typename spacing_func_t>
    void operator()(domain_t& domain, DomainShape<param_vec_t>& param_domain_shape,
                    const param_func_t& param_function, const jm_func_t& param_jacobian,
                    const spacing_func_t& spacing_function, int type = 0) const;

    /// Overload for constant function and Shape instead of DomainDiscretization
    template <typename param_func_t, typename jm_func_t>
    void operator()(domain_t& domain, DomainShape<param_vec_t>& param_domain_shape,
                    const param_func_t& param_function, const jm_func_t& param_jacobian,
                    const scalar_t& h, int type = 0) const {
        this->operator()(domain, param_domain_shape, param_function, param_jacobian,
                         [=] (const vec_t&) { return h; }, type);
    }
};
}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_GENERALSURFACEFILL_FWD_HPP_
