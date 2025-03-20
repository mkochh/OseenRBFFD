#ifndef MEDUSA_BITS_DOMAINS_GRIDFILL_HPP_
#define MEDUSA_BITS_DOMAINS_GRIDFILL_HPP_

/**
 * @file
 * Implementation of GridFill engine.
 */

#include "GridFill_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/spatial_search/KDTree.hpp>
#include <medusa/bits/utils/numutils.hpp>
#include <medusa/bits/domains/DomainDiscretization.hpp>

namespace mm {

template <typename vec_t>
GridFill<vec_t>::GridFill(const vec_t& bot, const vec_t& top) : bot(bot), top(top) {
    for (int i = 0; i < dim; ++i) {
        assert_msg(bot[i] <= top[i], "Bottom bounds %s must be below top bounds %s.", bot, top);
    }
}

template <typename vec_t>
void GridFill<vec_t>::operator()(DomainDiscretization<vec_t>& domain, const scalar_t& h,
                                 int type) const {
    if (type == 0) type = 1;
    Vec<int, dim> counts;
    for (int i = 0; i < dim; ++i) {
        counts[i] = iceil((top[i] - bot[i]) / h) + 1;
    }

    KDTree<vec_t> tree(domain.positions());
    for (const auto& x : linspace(bot, top, counts)) {
        if (domain.contains(x) && (tree.size() == 0 || tree.query(x, 1).second[0] >= h*h)) {
            domain.addInternalNode(x, type);
        }
    }
}

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_GRIDFILL_HPP_
