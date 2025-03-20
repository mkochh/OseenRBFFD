#ifndef MEDUSA_BITS_DOMAINS_HALFLINKSREFINE_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_HALFLINKSREFINE_FWD_HPP_

/**
 * @file
 * Declaration of the half-link refinement algorithm.
 *
 * @example test/domains/HalfLinksRefine_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/types/Range_fwd.hpp>
#include "DomainDiscretization_fwd.hpp"

namespace mm {

template <class vec_t>
class KDTreeMutable;

/**
 * Refine a region of nodes `region` by connecting every node in `region` to its support
 * domain and generating new nodes at half distances. The new nodes are filtered to meet a
 * minimum distance criterion to prevent points that would be too close.
 *
 * Usage example:
 * @snippet domains/HalfLinksRefine_test.cpp Half links usage example
 * @ingroup domains
 */
class HalfLinksRefine {
    Range<int> region_ = {};  ///< Node indexes around which to refine.
    double fraction_ = 0.4;  ///< Minimal distance fraction.

  public:
    /**
     * Set region to refine. The region is a set of indices of the nodes to be refined.
     * All nodes are refined by default.
     */
    HalfLinksRefine& region(Range<int> region) { region_ = std::move(region); return *this; }
    /**
     * Minimal distance criterion around point `p` is that nodes generated at `p` must be at least
     * `fraction * closest_support_node` away from all nodes. Fraction must be less than `1/2`.
     * Default value is HalfLinksRefine::fraction_;
     */
    HalfLinksRefine& min_dist(double fraction) { fraction_ = fraction; return *this; }

    /**
     * Refines given domain.
     * @return The indexes of the added nodes in `positions`.
     */
    template <typename vec_t>
    Range<int> operator()(DomainDiscretization<vec_t>& domain) const;

    /**
     * Refine the domain with already given tree.
     * @param domain Domain discretization object, which will be refined.
     * @param tree A mutable KDTree, which contains all domain nodes. New nodes will be inserted.
     * @return The indexes of the added nodes in `positions`.
     */
    template <typename vec_t>
    Range<int> operator()(DomainDiscretization<vec_t>& domain, KDTreeMutable<vec_t>& tree) const;

  private:
    /**
     * Refine implementation with more control over fine grained details.
     * @param domain Domain to refine.
     * @param region A list of indices to refine.
     * @param fraction Minimal distance fraction.
     * @param domain_tree Tree containing all domain points.
     * @return List of indices of added nodes.
     */
    template <class vec_t>
    static Range<int> refine_impl(
            DomainDiscretization<vec_t>& domain, const Range<int>& region,
            double fraction, KDTreeMutable<vec_t>& domain_tree);
};


}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_HALFLINKSREFINE_FWD_HPP_
