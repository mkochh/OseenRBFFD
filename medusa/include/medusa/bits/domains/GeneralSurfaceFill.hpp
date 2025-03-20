#ifndef MEDUSA_BITS_DOMAINS_GENERALSURFACEFILL_HPP_
#define MEDUSA_BITS_DOMAINS_GENERALSURFACEFILL_HPP_

/**
 * @file
 * Implementation of the general node placing algorithm for parametrically given surfaces.
 */

#include "GeneralSurfaceFill_fwd.hpp"
#include "discretization_helpers.hpp"
#include "compute_normal_fwd.hpp"

#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/types/Range.hpp>
#include <medusa/bits/utils/randutils.hpp>
#include <medusa/bits/spatial_search/KDTreeMutable.hpp>
#include <medusa/bits/spatial_search/KDTree.hpp>
#include <random>

namespace mm {

template <typename vec_t, typename param_vec_t>
GeneralSurfaceFill<vec_t, param_vec_t>::GeneralSurfaceFill() : seed_(get_seed()) {}

template <typename vec_t, typename param_vec_t>
template <typename param_func_t, typename jm_func_t, typename spacing_func_t>
void GeneralSurfaceFill<vec_t, param_vec_t>::operator()(domain_t& domain,
        param_domain_t& param_domain, const param_func_t& param_function,
        const jm_func_t& param_jacobian, const spacing_func_t& spacing_function, int type) const {
    assert_msg(domain.size() >= param_domain.size(),
               "Domain must contain all points corresponding to the parameters in the "
               "parametric domain. ");
    KDTreeMutable<vec_t> tree(domain.positions());
    operator()(domain, param_domain, param_function, param_jacobian, spacing_function, tree, type);
}
template <typename vec_t, typename param_vec_t>
template <typename param_func_t, typename jm_func_t, typename spacing_func_t>
void GeneralSurfaceFill<vec_t, param_vec_t>::operator()(domain_t& domain,
        DomainShape<param_vec_t>& param_domain_shape, const param_func_t& param_function,
        const jm_func_t& param_jacobian, const spacing_func_t& spacing_function, int type) const {
    DomainDiscretization<param_vec_t> param_domain(param_domain_shape);
    operator()(domain, param_domain, param_function, param_jacobian, spacing_function, type);
}

template <typename vec_t, typename param_vec_t>
template <typename param_func_t, typename jm_func_t, typename search_t, typename spacing_func_t>
void GeneralSurfaceFill<vec_t, param_vec_t>::operator()(domain_t& domain,
        param_domain_t& param_domain, const param_func_t& param_function,
        const jm_func_t& param_jacobian, const spacing_func_t& spacing_function, search_t& tree,
        int type) const {
    auto param = [&param_function, &param_jacobian](param_vec_t t) {
        return std::make_pair(param_function(t), param_jacobian(t));
    };
    fillParametrization(domain, param_domain, param, spacing_function, tree, type);
}

template <typename vec_t, typename param_vec_t>
template <typename param_func_t, typename search_t, typename spacing_func_t>
void GeneralSurfaceFill<vec_t, param_vec_t>::fillParametrization(domain_t& domain,
        param_domain_t& param_domain, const param_func_t& param,
        const spacing_func_t& spacing_function, search_t& tree, int type) const {
    if (type == 0) type = -1;
    std::mt19937 gen(seed_);
    KDTree<param_vec_t> param_boundary_search;
    param_domain.makeDiscreteContainsStructure(param_boundary_search);

    int cur_node = 0;
    int end_node = param_domain.size();
    if (end_node == 0) {
        assert_msg(param_domain.shape().hasContains(),
                "Parametric domain shape must have `contains` method implemented if empty "
                "parametric domain is given. Try supplying an initial point.");
        // If parameter domain is empty, pick a random node inside it.
        param_vec_t lo_bound, hi_bound, random_node;
        std::tie(lo_bound, hi_bound) = param_domain.shape().bbox();
        std::vector<std::uniform_real_distribution<scalar_t>> distributions;
        for (int j = 0; j < param_dim; ++j) distributions.emplace_back(lo_bound[j], hi_bound[j]);
        int count = 0;
        scalar_t d_sq, check_radius;
        vec_t point;
        Eigen::Matrix<scalar_t, dim, param_dim> jm;
        do {
            for (int j = 0; j < param_dim; ++j) { random_node[j] = distributions[j](gen); }
            // Node must be at least spacing_function(random_node) away from other nodes
            std::tie(point, jm) = param(random_node);
            check_radius = spacing_function(point);
            d_sq = tree.size() == 0 ? 10 * check_radius * check_radius :
                                      tree.query(point).second[0];
            if (++count >= 10000) {
                std::string message = "No suitable node in parametric domain could be found after "
                                      "10000 tries. This might happen if parametric domain volume "
                                      "is very small compared to the volume of the bounding box or"
                                      " if there are a lot of nodes in input domain.";
                throw std::runtime_error(message);
            }
        } while (!param_domain.contains(random_node, param_boundary_search) ||
                    !(d_sq >= (zeta * check_radius) * (zeta * check_radius)));
        param_domain.addInternalNode(random_node, 1);
        vec_t normal = surface_fill_internal::compute_normal(jm);
        tree.insert(point);
        domain.addBoundaryNode(point, type, normal);
        end_node = 1;
    }

    // Main algorithm loop.
    while (cur_node < end_node && end_node < max_points) {
        param_vec_t param_point = param_domain.pos(cur_node);
        vec_t initial_pt;
        Eigen::Matrix<scalar_t, dim, param_dim> jm;
        std::tie(initial_pt, jm) = param(param_point);

        auto unit_candidates = discretization_helpers::SphereDiscretization<scalar_t, param_dim>::
                               construct(1, n_samples, gen);

        // Filter unit_candidates regarding the domain and proximity criteria.
        for (const auto& u_cand : unit_candidates) {
            scalar_t alpha = spacing_function(initial_pt) / (jm * u_cand).norm();
            param_vec_t node = param_point + alpha * u_cand;  // Shift center to point `param_point`
            if (!param_domain.contains(node, param_boundary_search)) continue;

            vec_t new_pt;
            Eigen::Matrix<scalar_t, dim, param_dim> new_jm;
            std::tie(new_pt, new_jm) = param(node);

            scalar_t d_sq = tree.query(new_pt).second[0];
            // Check radius must be new radius, otherwise algorithm might terminate in 2d.
            scalar_t check_radius = (new_pt - initial_pt).norm();
            // Check if node is far enough from other nodes.
            if (d_sq < (zeta * check_radius) * (zeta * check_radius)) continue;

            // Insert into param domain.
            param_domain.addInternalNode(node, 1);
            // Insert into kd tree.
            tree.insert(new_pt);
            // Inser into domain.
            vec_t normal = surface_fill_internal::compute_normal(jm);
            domain.addBoundaryNode(new_pt, type, normal);

            end_node++;
        }
        cur_node++;
    }
    if (end_node >= max_points) {
        std::cerr << "Maximum number of points reached, fill may be incomplete." << std::endl;
    }
}

template <typename vec_t, typename param_vec_t>
GeneralSurfaceFill<vec_t, param_vec_t>& GeneralSurfaceFill<vec_t, param_vec_t>::proximityTolerance(
        scalar_t zeta) {
    assert_msg((0 < zeta && zeta < 1), "Zeta must be between 0 and 1, got %f.", zeta);
    this->zeta = zeta;
    return *this;
}
}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_GENERALSURFACEFILL_HPP_
