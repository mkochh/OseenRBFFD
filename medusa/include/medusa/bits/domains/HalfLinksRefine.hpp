#ifndef MEDUSA_BITS_DOMAINS_HALFLINKSREFINE_HPP_
#define MEDUSA_BITS_DOMAINS_HALFLINKSREFINE_HPP_

/**
 * @file
 * Implementation of the half-link refinement algorithm.
 */


#include "HalfLinksRefine_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/spatial_search/KDTreeMutable.hpp>
#include <medusa/bits/types/Range.hpp>
#include <numeric>

namespace mm {

template <typename vec_t>
Range<int> HalfLinksRefine::operator()(DomainDiscretization<vec_t>& domain) const {
    KDTreeMutable<vec_t> domain_tree(domain.positions());
    return operator()(domain, domain_tree);
}

template <typename vec_t>
Range<int> HalfLinksRefine::operator()(
        DomainDiscretization<vec_t>& domain, KDTreeMutable<vec_t>& tree) const {
    Range<int> region = region_;
    if (region.empty()) { region = domain.all(); }

    // sort: boundary nodes first
    std::sort(region.begin(), region.end(),
              [&](int i, int j) { return domain.type(i) < domain.type(j); });

    return refine_impl(domain, region, fraction_, tree);
}

template <class vec_t>
Range<int> HalfLinksRefine::refine_impl(
        DomainDiscretization<vec_t>& domain, const Range<int>& region,
        double fraction, KDTreeMutable<vec_t>& domain_tree) {
    int region_size = region.size();
    int N = domain.size();
    using scalar_t = typename vec_t::scalar_t;

    assert_msg(region_size > 0, "The region to refine is empty.");

    // Iterate through points in region and generate new points
    int num_new_points = 0;
    for (int i = 0; i < region_size; i++) {
        int c = region[i];  // the global domain index
        const vec_t pos = domain.pos(c);
        const Range<int> supp = domain.support(c);
        assert_msg(supp.size() >= 2, "At least 2 nodes must be in support of every node, %d "
                                     "found in support of node %d.", supp.size(), c);
        scalar_t min_dist = fraction * domain.dr(c);

        int n = supp.size();
        // Half links to my support.
        for (int j = 1; j < n; ++j) {
            int s = supp[j];
            vec_t candidate = 0.5 * (domain.pos(c) + domain.pos(s));

            // Decide the type of new node and project to boundary if necessary.
            if (domain.type(c) == 0 || domain.type(s) == 0) continue;
            int candidate_type = std::max(domain.type(c), domain.type(s));
            vec_t normal_vector;
            if (candidate_type < 0) {
                normal_vector = domain.normal(c) + domain.normal(s);
                if (normal_vector.squaredNorm() < 1e-15) {
                    // normal_vector of given points point in opposite directions.
                    continue;
                }
                normal_vector.normalize();
                auto result = domain.shape().projectPointToBoundary(candidate, normal_vector);
                if (result.first) {
                    candidate = result.second;
                } else {
                    std::cerr << format("Adding point %s with type %d to the boundary along "
                                        "normal %s was not successful.",
                                        candidate, candidate_type, normal_vector)
                              << std::endl;
                    continue;
                }
            }

            // If nodes are out of the domain, they should be thrown away.
            if (!domain.shape().contains(candidate)) continue;
            // Distance to closest node.
            double dist2 = domain_tree.query(candidate).second[0];
            if (dist2 >= min_dist * min_dist) {  // add new point
                ++num_new_points;
                domain_tree.insert(candidate);
                if (candidate_type < 0) {
                    domain.addBoundaryNode(candidate, candidate_type, normal_vector);
                } else {
                    domain.addInternalNode(candidate, candidate_type);
                }
            }
        }
    }

    // return back indices of new nodes
    Range<int> new_ids(num_new_points);
    std::iota(new_ids.begin(), new_ids.end(), N);
    return new_ids;
}

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_HALFLINKSREFINE_HPP_
