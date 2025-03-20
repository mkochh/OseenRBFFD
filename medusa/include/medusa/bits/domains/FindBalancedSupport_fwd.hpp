#ifndef MEDUSA_BITS_DOMAINS_FINDBALANCEDSUPPORT_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_FINDBALANCEDSUPPORT_FWD_HPP_

/**
 * @file
 * Declaration of class for finding support nodes.
 *
 * @example test/domains/FindBalancedSupport_test.cpp
 */

#include "Eigen/Core"
#include <medusa/Config.hpp>
#include <medusa/bits/types/Range_fwd.hpp>

namespace mm {

/**
 * Class representing the engine for finding directionally balanced supports.
 * The support is directionally balanced, if it has nodes in all axis aligned
 * hyper-quadrants centered around a node. It the node is a boundary node,
 * then the same must be satisfied for half space, defined by the normal.
 *
 * Usage example:
 * @snippet domains/FindBalancedSupport_test.cpp FindBalancedSupport usage example
 * @ingroup domains
 */
class FindBalancedSupport {
  public:
    int min_support_;  ///< Minimal support size.
    int max_support_;  ///< Maximal support size.
    Range<int> for_which_;  ///< Find support only for these nodes.
    Range<int> search_among_;  ///< Search only among these nodes.
    bool force_self_;  ///< Force each node as the first element of its support.

  public:
    /// Constructs an engine with given min and max support sizes.
    FindBalancedSupport(int min_support, int max_support);
    /**
     * Find support only for these nodes. If not given, finds support for all non-zero type
     * domain nodes, as defined by @ref DomainDiscretization::all "domain.all()".
     */
    FindBalancedSupport& forNodes(indexes_t for_which);
    /**
     * Search only among given nodes. If not given, searches among all nodes, including
     * zero-type nodes.
     */
    FindBalancedSupport& searchAmong(indexes_t search_among);
    /// Put each node as the first of its support, even if it is not included in searchAmong().
    FindBalancedSupport& forceSelf(bool b = true);
    /// Set minimum support size. This overrides the size set in constructor.
    FindBalancedSupport& minSupportSize(int size);
    /// Set maximal support size. This overrides the size set in constructor.
    FindBalancedSupport& maxSupportSize(int size);

    /// Find support for nodes in domain.
    template <typename domain_t>
    void operator()(domain_t& domain) const;

    /// Finds support for a given node `i`.
    template <typename domain_t, typename kdtree_t>
    static Range<int> balancedSupport(
            domain_t& domain, const kdtree_t& tree, int i, int min_support_size,
            int max_support_size, const Range<int>& search_among, bool force_self);

  private:
    /// Gets orthonormal basis of vectors perpendicular to `normal`.
    template <class vec_t>
    static Eigen::Matrix<typename vec_t::scalar_t, vec_t::dim, vec_t::dim-1>
    getFrame(const vec_t& normal);

    /**
     * Marks a quadrant in which dx resides as occupied.
     * Returns `true` if this quadrant is newly occupied and `false` otherwise.
     */
    template <class vec_t>
    static bool mark_quadrant(const vec_t& dx, typename vec_t::scalar_t tol,
                              std::vector<bool>& octants_covered,
                              const Eigen::Matrix<typename vec_t::scalar_t, vec_t::dim, -1>& basis);
};

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_FINDBALANCEDSUPPORT_FWD_HPP_
