#ifndef MEDUSA_BITS_DOMAINS_BASICRELAX_HPP_
#define MEDUSA_BITS_DOMAINS_BASICRELAX_HPP_

#include "BasicRelax_fwd.hpp"
#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/types/Range.hpp>
#include <medusa/bits/spatial_search/KDTree.hpp>
#include <medusa/bits/spatial_search/KDTreeMutable.hpp>

/**
 * @file
 * Implementation of class for box shaped domains.
 */

namespace mm {


template <class domain_t, class radius_func_type>
void BasicRelax::operator()(domain_t& domain, const radius_func_type& r_func) const {
    typedef typename domain_t::vector_t vec_t;
    typedef typename domain_t::scalar_t scalar_t;
    const int dim = domain_t::dim;

    // if no nodes are supplied, all interior nodes are processed
    Range<int> nodes = (nodes_.empty()) ? static_cast<Range<int>>(domain.interior()) : nodes_;

    int N = domain.size();
    assert_msg(N > num_neighbours, "Cannot relax a domain with less that num_neighbours nodes, "
                                   "num_neighbours = %d, N = %d.", num_neighbours, N);
    assert_msg(num_neighbours >= 1, "At least two neighbours are required, got %d.",
               num_neighbours);
    for (int i : nodes) {
        assert_msg(0 <= i && i < N, "Node index %d out of range [0, %d) in relax.", i, N);
    }

    Range<int> bnd = domain.boundary();
    assert_msg(!bnd.empty(), "Relax requires boundary nodes, but this domain has none!");
    KDTreeMutable<vec_t> boundary_tree(domain.positions()[bnd]);
    KDTree<vec_t> tree;
    Range<int> indices;

    Range<int> to_delete = {};  // nodes to be deleted from domain after relax
    bool removed_nodes_warning = false;

    for (int c = 0; c < num_iterations; ++c) {
        if (c == 0 || c % rebuild_tree_after == 0) {
            tree.reset(domain.positions());
        }
        // nodes to be removed from relax procedure
        Range<int> to_remove = {};
        // Annealing factor -- depends only on iteration step
        double a_f = ((initial_heat - final_heat) * (num_iterations - c) / num_iterations
                      + final_heat);
        // Range of all displacements
        Range<vec_t> dp(nodes.size(), vec_t(0.0));
        int k;
//            #if defined(_OPENMP)
//            #pragma omp parallel for private(k) schedule(static)
        // before enabling openmp vectors to_delete and to_remove
        // have to be introduced for each thread and combined at the end of the loop
        // otherwise expect se.g.\ fault
//            #endif
        for (k = 0; k < nodes.size(); ++k) {
            int i = nodes[k];
            assert_msg(domain.type(i) > 0, "Only interior nodes are to be relaxed");
            vec_t pos = domain.pos(i);
            indices = tree.query(pos, num_neighbours + 1).first;
            double d0 = std::pow((domain.pos(indices[1]) - pos).norm(), 1.0 / dim);
            // loop over support and compute potential contributions
            scalar_t r_0 = 0;
            for (int j = 1; j < indices.size(); ++j) {
                // central node to support nodes vector
                vec_t r = domain.pos(indices[j]) - pos;
                assert_msg(r.norm() > 1e-15, "Nodes %d and %d have the same coordinates %s. "
                                             "Try using lower initial heat.", i, indices[j], pos);
                // potential gradient
                // target density in terms of radius to the closest node
                r_0 = r_func(domain.pos(indices[j]));
                dp[k] += std::pow(std::abs(r_0) / r.norm(), potential_order) * r;
            }
            // no-distribution mode
            if (r_0 < 0) dp[k] = 0.2 * dp[k].normalized() * a_f * d0;
            dp[k] = -a_f * d0 * dp[k];
            // Project escaped nodes on the boundary
            if (!domain.shape().contains(pos + dp[k])) {
                vec_t tmp, normal;  // projection temp variables
                bool project = true;
                switch (projection_type) {
                    case DO_NOT_PROJECT:  // mark escaped nodes for deletion
                        to_remove.push_back(k);
                        to_delete.push_back(nodes[k]);
                        project = false;
                        break;
                    case PROJECT_IN_DIRECTION:  // project on boundary in direction of relax
                        tmp = 0.5 * (2 * pos + dp[k]);
                        normal = dp[k].normalized();
                        break;
                    case PROJECT_BETWEEN_CLOSEST:  // project between closest nodes on boundary
                        Range<int> idx = boundary_tree.query(pos, 2).first;
                        int s = bnd[idx[0]], t = bnd[idx[1]];
                        tmp = 0.5 * (domain.pos(s) + domain.pos(t));
                        normal = (0.5 * (domain.normal(s) + domain.normal(t))).normalized();
                        break;
                }  // projection method selection end
                if (project) {
                    // performs actual projection
                    bool success;
                    vec_t proj;
                    std::tie(success, proj) = domain.shape().projectPointToBoundary(tmp, normal);
                    if (success) {
                        // check if projected node is not too close to boundary node
                        Range<double> dist;
                        Range<int> ind;
                        std::tie(ind, dist) = boundary_tree.query(proj, 2);
                        if (std::sqrt(dist[0] / dist[1]) < boundary_projection_threshold) {
                            dp[k] = -dp[k];
                            continue;
                        }
                        double d_1 = std::sqrt(dist[0]);
                        double d_2 = std::sqrt(dist[1]);
                        double d_0 = std::sqrt(dist[0] + dist[1]);

                        normal = ((1 - d_1 / d_0) * domain.normal(bnd[ind[0]])
                                  + (1 - d_2 / d_0) * domain.normal(bnd[ind[1]])).normalized();
                        int boundary_type = domain.type(bnd[ind[0]]);
                        domain.changeToBoundary(i, proj, boundary_type, normal);
                        boundary_tree.insert(proj);
                        bnd.push_back(i);
                        // mark projected node to be removed from the list for relax
                        to_remove.push_back(k);
                    } else {
                        removed_nodes_warning = true;
                        // if projection failed, remove node
                        to_remove.push_back(k);
                        to_delete.push_back(nodes[k]);
                    }
                }  // project if end
            }  // escaped nodes if end
        }  // spatial loop end
        nodes.remove(to_remove);
        assert_msg(nodes.size() > num_neighbours + 1,
                   "No nodes in relax pool anymore, perhaps use lower heat");
        // apply displacement
        // #if defined(_OPENMP)
        //    #pragma omp parallel for private(k) schedule(static)
        // #endif
        for (k = 0; k < nodes.size(); ++k) {
            domain.pos(nodes[k]) += dp[k];
        }
        // DEBUG OUT -- for plotting
        // debug_out << "nodes(" << c + 1 << ",:,:)=" << domain.positions << ";\n";
    }
    // remove nodes at the end of relax to avoid mess with indices
    domain.removeNodes(to_delete);
    #if 0  // debug
    debug_out << "nodes_final(:,:)=" << domain.positions << ";\n";
        debug_out << "normal =" << domain.normals() << ";\n";
        Range<int> boundary = domain.types < 0;
        debug_out << "boundary =" << boundary << ";\n";
        debug_out << "boundary_map =" << domain.boundary_map() << ";\n";
        debug_out.close();
    #endif
    if (removed_nodes_warning) {
        print_red("warning -- nodes removed during relax due to projection fail\n");
    }
}

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_BASICRELAX_HPP_
