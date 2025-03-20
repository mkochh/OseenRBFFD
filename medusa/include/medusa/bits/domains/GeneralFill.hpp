#ifndef MEDUSA_BITS_DOMAINS_GENERALFILL_HPP_
#define MEDUSA_BITS_DOMAINS_GENERALFILL_HPP_

/**
 * @file
 * Implementation of the general node placing algorithm.
 */

#include "GeneralFill_fwd.hpp"
#include "discretization_helpers.hpp"

#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/types/Range.hpp>
#include <random>
#include <medusa/bits/utils/randutils.hpp>
#include <medusa/bits/spatial_search/KDTreeMutable.hpp>
#include <medusa/bits/spatial_search/KDTree.hpp>

namespace mm {

template <typename vec_t>
GeneralFill<vec_t>::GeneralFill() : seed_(get_seed()) {}

/// @cond
template <typename vec_t>
template <typename func_t>
void GeneralFill<vec_t>::operator()(DomainDiscretization<vec_t>& domain, const func_t& h,
                                    int type) const {
    KDTreeMutable<vec_t> tree;
    operator()(domain, h, tree, type);
}

template <typename vec_t>
template <typename func_t, typename search_structure_t>
void GeneralFill<vec_t>::operator()(DomainDiscretization<vec_t>& domain, const func_t& h,
                                    search_structure_t& search, int type) const {
    KDTree<vec_t> contains_search;
    domain.makeDiscreteContainsStructure(contains_search);
    auto contains_function = [&] (const vec_t point) {
        return domain.contains(point, contains_search);
    };
    operator()(domain, h, search, contains_function, type);
}

template <typename vec_t>
template <typename func_t, typename search_structure_t, typename contains_function_t>
void GeneralFill<vec_t>::operator()(DomainDiscretization<vec_t>& domain, const func_t& h,
                                    search_structure_t& search,
                                    contains_function_t& contains_function, int type) const {
    if (type == 0) type = 1;
    std::mt19937 gen(seed_);

    int cur_node = 0;
    int end_node = domain.size();
    if (end_node == 0) {  // If domain is empty, pick a random node inside it.
        vec_t lo_bound, hi_bound, random_node;
        std::tie(lo_bound, hi_bound) = domain.shape().bbox();
        std::vector<std::uniform_real_distribution<scalar_t>> distributions;
        for (int j = 0; j < dim; ++j) distributions.emplace_back(lo_bound[j], hi_bound[j]);
        int count = 0;
        do {
            for (int j = 0; j < dim; ++j) { random_node[j] = distributions[j](gen); }
            if (++count >= 10000) {
                std::string message = "No suitable node in domain could be found after 10000 tries."
                                      " This might happen if domain volume is very small compared "
                                       "to the volume of the bounding box. Try manually supplying "
                                      "the initial point.";
                throw std::runtime_error(message);
            }
        } while (!contains_function(random_node));
        domain.addInternalNode(random_node, type);
        end_node = 1;
    }

    // Main algorithm loop.
    search.insert(domain.positions());
    while (cur_node < end_node && end_node < max_points) {
        vec_t p = domain.pos(cur_node);
        scalar_t r = h(p);
        assert_msg(r > 0, "Nodal spacing radius should be > 0, got %g.", r);

        auto candidates = discretization_helpers::SphereDiscretization<scalar_t, dim>
                                                ::construct(r, n_samples, gen);
        // filter candidates regarding the domain and proximity criteria
        for (const auto& f : candidates) {
            vec_t node = p + f;  // Shift center to point `p`.
            if (!contains_function(node)) continue;  // If node is not in the domain.
            if (search.existsPointInSphere(node, r*zeta)) continue;  // If node is not far enough
            domain.addInternalNode(node, type);                      // from other nodes.
            search.insert(node);
            end_node++;
        }
        cur_node++;
    }
    if (end_node >= max_points) {
        std::cerr << "Maximum number of points reached, fill may be incomplete." << std::endl;
    }
}

template <typename vec_t>
GeneralFill<vec_t>& GeneralFill<vec_t>::proximityTolerance(scalar_t zeta) {
    assert_msg((0 < zeta && zeta < 1), "Zeta must be between 0 and 1, got %f.", zeta);
    this->zeta = zeta;
    return *this;
}
/// @endcond

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_GENERALFILL_HPP_
