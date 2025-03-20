#ifndef MEDUSA_BITS_DOMAINS_FINDCLOSEST_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_FINDCLOSEST_FWD_HPP_

/**
 * @file
 * Declarations of FindClosest.
 *
 * @example test/domains/FindClosest_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/types/Range_fwd.hpp>

namespace mm {

/**
 * Class representing the engine for finding supports consisting of closest nodes.
 *
 * Usage example:
 * @snippet domains/FindClosest_test.cpp FindClosest usage example
 * @ingroup domains
 */
class FindClosest {
  private:
    int support_size;  ///< Support size.
    Range<int> for_which_;  ///< Find support only for these nodes.
    Range<int> search_among_;  ///< Search only among these nodes.
    bool force_self_;  ///< Force each node as the first element of its support.

  public:
    /// Constructs an engine for finding support with given support size.
    FindClosest(int support_size) : support_size(support_size), for_which_(), search_among_(),
            force_self_(false) {}
    /**
     * Find support only for these nodes. If not given, finds support for all non-zero type
     * domain nodes, as defined by @ref DomainDiscretization::all "domain.all()".
     */
    FindClosest& forNodes(indexes_t for_which) { for_which_ = std::move(for_which); return *this; }
    /**
     * Search only among given nodes. If not given, searches among all nodes, including
     * zero-type nodes.
     */
    FindClosest& searchAmong(indexes_t search_among) {
        search_among_ = std::move(search_among); return *this; }
    /// Put each node as the first of its support, even if it is not included in searchAmong().
    FindClosest& forceSelf(bool b = true) { force_self_ = b; return *this; }
    /// Find `num` closest nodes. This methods overrides the value set in the constructor.
    FindClosest& numClosest(int num) { support_size = num; return *this; }

    /// Find support for nodes in domain.
    template <typename domain_t>
    void operator()(domain_t& domain) const;
};

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_FINDCLOSEST_FWD_HPP_
