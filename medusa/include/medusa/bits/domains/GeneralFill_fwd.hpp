#ifndef MEDUSA_BITS_DOMAINS_GENERALFILL_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_GENERALFILL_FWD_HPP_

/**
 * @file
 * Declaration of the general node placing algorithm.
 *
 * @example test/domains/GeneralFill_test.cpp
 */

#include <medusa/Config.hpp>

/**
 * @file
 * Declaration of GeneralFill class.
 */

namespace mm {

template <typename vec_t>
class DomainDiscretization;

/**
 * Implements general `n`-d node placing algorithm, as described in https://arxiv.org/abs/1812.03160
 * If you specifically use this algorithm, we would appreciate if you cite the above publication.
 *
 * The algorithm starts from existing nodes, usually the boundary discretization nodes.
 * They are processed in a queue and for each node, new candidates are generated around it on
 * a circle with radius defined by a supplied density function. Candidates are then checked
 * and some are accepted and added to the queue. This is repeated until no more nodes
 * can be generated or the maximal number of points has been reached.
 *
 * Usage example:
 * @snippet domains/GeneralFill_test.cpp PDS usage example
 *
 * Using background grid:
 * @snippet domains/GeneralFill_test.cpp Background grid
 *
 * @sa GrainDropFill, BasicRelax
 * @ingroup domains
 */
template <typename vec_t>
class GeneralFill {
  public:
    typedef typename vec_t::scalar_t scalar_t;  ///< Scalar type;
    typedef vec_t vector_t;  ///< Vector type.
    /// Store dimension of the domain.
    enum { /** Dimensionality of the domain. */ dim = vec_t::dim };

  private:
    int max_points = 5000000;  ///< Maximal number of points generated.
    int seed_;  ///<  Seed for the random number generator.
    int n_samples = 15;  ///< Number of samples.
    scalar_t zeta = 1 - 1e-10;  ///< Proximity tolerance.

  public:
    GeneralFill();

    /** Maximal number of points generated
     * @warning The actual upper bound for number of points generated is maxPoints @f$\pm@f$
     * numSamples.
     */
    GeneralFill& maxPoints(int max_points) { this->max_points = max_points; return *this; }
    /// Set custom seed for the random number generator.
    GeneralFill& seed(int seed) { seed_ = seed; return *this; }
    /// Set proximity tolerance. A new candidate mush be at least `zeta*h(p)` away. It should hold
    /// that `0 < zeta < 1`.
    GeneralFill& proximityTolerance(scalar_t zeta);
    /**
     * Controls the number of generated candidates from each point. For 2-D it is the actual number
     * of candidates, for 3-D it is the number of candidates on the great circle. Its value
     * is ignored in 1-D.
     */
    GeneralFill& numSamples(int n_samples) { this->n_samples = n_samples; return *this; }

    /// Fills given domain according to the nodal spacing function `h`.
    template <typename func_t>
    void operator()(DomainDiscretization<vec_t>& domain, const func_t& h, int type = 0) const;

    /**
     * Fills domain with a quality node distribution.
     * @param domain Domain to fill with nodes. Can be partially filled.
     * @param h Nodal spacing function.
     * @param search Spatial search structure to be used, e.g. a KDTree or KDGrid.
     * @param contains_function Function taking a point and returning `True` if it is contained in
     * the domain and `False` otherwise.
     * @param type Type of the nodes. If 0 the engines default value is used.
     *
     * @warning Header <code> \#include <medusa/bits/domains/GeneralFill.hpp></code> must be
     * included additionally to <code> \#include <medusa/Medusa_fwd.hpp></code>.
     */
    template <typename func_t, typename search_structure_t, typename contains_function_t>
    void operator()(DomainDiscretization<vec_t>& domain, const func_t& h,
                    search_structure_t& search, contains_function_t& contains_function,
                    int type = 0) const;


    /// Overload for constant function.
    template <typename search_structure_t, typename contains_function_t>
    void operator()(DomainDiscretization<vec_t>& domain, const scalar_t& h,
            search_structure_t& search, contains_function_t& contains_function,
            int type = 0) const {
        this->operator()(domain, [=] (const vector_t&) { return h; }, search, contains_function,
                type);
    }

    /// Overload for constant function and default search and contains structures.
    void operator()(DomainDiscretization<vec_t>& domain, const scalar_t& h, int type = 0) const {
        this->operator()(domain, [=] (const vector_t&) { return h; }, type);
    }

    /// Overload without providing contains_function, using the default
    /// DomainDiscretization::contains.
    template <typename func_t, typename search_structure_t>
    void operator()(DomainDiscretization<vec_t>& domain, const func_t& h,
                    search_structure_t& search, int type = 0) const;
};

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_GENERALFILL_FWD_HPP_
