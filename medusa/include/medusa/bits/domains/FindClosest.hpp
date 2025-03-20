#ifndef MEDUSA_BITS_DOMAINS_FINDCLOSEST_HPP_
#define MEDUSA_BITS_DOMAINS_FINDCLOSEST_HPP_

/**
 * @file
 * Implementation of FindClosest.
 */

#include "FindClosest_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/spatial_search/KDTree.hpp>

/**
 * @file
 * Implementation of FindClosest class.
 */

namespace mm {

template <typename domain_t>
void FindClosest::operator()(domain_t& domain) const {
    // TODO(jureslak): speed up if search among is empty
    auto for_which = for_which_;
    if (for_which.empty()) for_which = domain.all();
    auto search_among = search_among_;
    if (search_among.empty()) {
        search_among = Range<int>::seq(domain.size());  // include potential ghost nodes
    }

    assert_msg(!domain.positions().empty(), "Cannot find support in an empty domain.");
    assert_msg(support_size > 0, "Support size must be greater than 0, got %d.", support_size);
    assert_msg(support_size <= search_among.size(),
               "Support size (%d) cannot exceed number of points that we are searching among (%d).",
               support_size, search_among.size());
    assert_msg(!for_which.empty(), "Set of nodes for which to find the support is empty. "
                                   "There appear to be no nodes in domain discretization.");
    assert_msg(!search_among.empty(), "Set of nodes to search among is empty. "
                                      "There appear to be no nodes in domain discretization.");
    for (int x : for_which) {
        assert_msg(0 <= x && x < domain.size(), "Index %d out of range [%d, %d) in forNodes.",
                   x, 0, domain.size());
    }
    for (int x : search_among) {
        assert_msg(0 <= x && x < domain.size(), "Index %d out of range [%d, %d) in "
                                                "searchAmong.", x, 0, domain.size());
    }
    KDTree<typename domain_t::vector_t> tree(domain.positions()[search_among]);
    if (force_self_) {
        for (int i : for_which) {
            Range<int> ss = search_among[tree.query(domain.pos(i), support_size).first];
            if (i != ss[0]) {
                domain.support(i) = {i};
                domain.support(i).insert(domain.support(i).end(), ss.begin(), ss.end()-1);
            } else {
                domain.support(i) = ss;
            }
        }
    } else {
        for (int i : for_which) {
            domain.support(i) = search_among[tree.query(domain.pos(i), support_size).first];
        }
    }
}

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_FINDCLOSEST_HPP_
