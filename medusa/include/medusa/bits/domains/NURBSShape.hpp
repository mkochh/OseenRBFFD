#ifndef MEDUSA_BITS_DOMAINS_NURBSSHAPE_HPP_
#define MEDUSA_BITS_DOMAINS_NURBSSHAPE_HPP_

/**
 * @file
 * Implementation of NURBSShape class.
 */

#include "NURBSShape_fwd.hpp"
#include "GeneralSurfaceFill.hpp"

#include <medusa/bits/domains/DomainDiscretization.hpp>
#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/spatial_search/KDTreeMutable.hpp>
#include <medusa/bits/domains/NURBSPatch.hpp>
#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/domains/compute_normal.hpp>

namespace mm {

template <typename vec_t, typename param_vec_t>
NURBSShape<vec_t, param_vec_t>::NURBSShape(const Range<NURBSPatch<vec_t, param_vec_t>>& patches_in)
                   : patches { patches_in }, seed_(get_seed()) {
    if (param_dim == 1) {
        n_samples = 2;
    }
    for (NURBSPatch<vec_t, param_vec_t>& patch : patches) {
        patch.computeDerivativeStructure();
    }
}

template <typename vec_t, typename param_vec_t>
NURBSShape<vec_t, param_vec_t>::NURBSShape(Range<NURBSPatch<vec_t, param_vec_t>>&& patches_in)
            : patches { patches_in }, seed_(get_seed()) {
    if (param_dim == 1) {
        n_samples = 2;
    }
    for (NURBSPatch<vec_t, param_vec_t> &patch : patches) {
        patch.computeDerivativeStructure();
    }
}

template <typename vec_t, typename param_vec_t>
DomainDiscretization<vec_t> NURBSShape<vec_t, param_vec_t>::discretizeBoundaryWithDensity(
        const std::function<scalar_t(vec_t)>& h, int type) const {
    return nurbs_shape_internal::NURBSShapeHelper<vec_t, param_vec_t>::
                                 discretizeBoundaryWithDensity(*this, h, type);
}

template <typename vec_t, typename param_vec_t>
NURBSShape<vec_t, param_vec_t>& NURBSShape<vec_t, param_vec_t>::proximityTolerance(
        scalar_t zeta) {
    assert_msg((0 < zeta && zeta < 1), "Zeta must be between 0 and 1, got %f.", zeta);
    this->zeta = zeta;
    return *this;
}

template <typename vec_t, typename param_vec_t>
NURBSShape<vec_t, param_vec_t>& NURBSShape<vec_t, param_vec_t>::boundaryProximity(
                                                                            scalar_t epsilon) {
    assert_msg((0 < epsilon && epsilon < 1), "Epsilon must be between 0 and 1, got %f.", epsilon);
    this->epsilon = epsilon;
    return *this;
}

/// @cond
template <typename vec_t>
struct nurbs_shape_internal::NURBSShapeHelper<vec_t, Vec<typename vec_t::scalar_t, 2>>{
    typedef Vec<typename vec_t::scalar_t, 2> param_vec_t;
    typedef typename vec_t::scalar_t scalar_t;
    enum { dim = vec_t::dim,
        proj_dim = dim + 1,
        param_dim = param_vec_t::dim};
    typedef Vec<scalar_t, proj_dim> proj_vec_t;

    static DomainDiscretization<vec_t> discretizeBoundaryWithDensity(
            const NURBSShape<vec_t, param_vec_t>& shape,
            const std::function<scalar_t(vec_t)>& h, int type) {
        if (type == 0) type = -1;
        std::mt19937 gen(shape.seed_);

        GeneralSurfaceFill<vec_t, param_vec_t> gsf;
        gsf.seed(shape.seed_).maxPoints(shape.max_points).proximityTolerance(shape.zeta)
           .numSamples(shape.n_samples);

        // Global domain and KDTree.
        DomainDiscretization<vec_t> domain(shape);
        KDTreeMutable<vec_t> tree_global;

        // Discretize every patch.
        for (const NURBSPatch<vec_t, param_vec_t>& patch : shape.patches) {
            BoxShape<param_vec_t> bs = patch.getDomain();
            DomainDiscretization<param_vec_t> local_param_d(bs);

            // Discretize boundaries first.
            Range<NURBSPatch<vec_t, Vec<scalar_t, 1>>> boundaries = patch.getBoundaries();
            int i = -1;
            for (NURBSPatch<vec_t, Vec<scalar_t, 1>>& boundary : boundaries) {
                i++;
                boundary.computeDerivativeStructure();

                KDTreeMutable<vec_t> tree_local;  // Local KD-Tree should be per-boundary.
                // Generate seed node.
                BoxShape<Vec<scalar_t, 1>> boundary_bs = boundary.getDomain();
                Vec<scalar_t, 1> t_min = boundary_bs.beg();
                Vec<scalar_t, 1> t_max = boundary_bs.end();

                std::uniform_real_distribution<> dis(t_min(0), t_max(0));
                Vec<scalar_t, 1> t_seed = Vec<scalar_t, 1>(dis(gen));

                // Insert into local structures.
                param_vec_t param_seed = patch.getPatchParameterFromBoundaryParameter(t_seed,
                                                                               i, shape.epsilon);
                scalar_t w_seed;
                vec_t pt_seed = boundary.evaluate(t_seed, &w_seed);
                local_param_d.addInternalNode(param_seed, 1);
                tree_local.insert(pt_seed);

                // Insert into global structures.
                scalar_t check_radius_seed = h(pt_seed);
                scalar_t d_sq = tree_global.size() == 0 ? 10 * check_radius_seed * check_radius_seed
                        : tree_global.query(pt_seed).second[0];
                if (d_sq >= (shape.zeta * check_radius_seed) * (shape.zeta * check_radius_seed)) {
                    vec_t normal_seed = surface_fill_internal::
                        compute_normal(patch.jacobian(param_seed, pt_seed, w_seed));
                    domain.addBoundaryNode(pt_seed, type, normal_seed);
                    tree_global.insert(pt_seed);
                }

                // General surface fill the boundary.
                int cur_node = local_param_d.size() - 1;
                int end_node = local_param_d.size();
                while (cur_node < end_node && end_node < shape.max_points) {
                    param_vec_t param = local_param_d.pos(cur_node);
                    Vec<scalar_t, 1> t = patch.getBoundaryParameterFromPatchParameter(param, i);
                    vec_t pt, der;
                    std::tie(pt, der) = boundary.evaluatePointAndJacobian(t);

                    std::array<Vec<scalar_t, 1>, 2> candidates{Vec<scalar_t, 1>(1),
                                                               Vec<scalar_t, 1>(-1)};

                    // Filter candidates.
                    scalar_t alpha = h(pt) / der.norm();
                    for (const auto& u_cand : candidates) {
                        Vec<scalar_t, 1> t_new = t + alpha * u_cand;  // Generate candidate.
                        param_vec_t param_new = patch.getPatchParameterFromBoundaryParameter(
                                t_new, i, shape.epsilon);
                        if (!local_param_d.contains(param_new)) continue;

                        scalar_t w_new;
                        vec_t pt_new = boundary.evaluate(t_new, &w_new);

                        // Check radius must be new radius.
                        scalar_t check_radius = (pt - pt_new).norm();

                        // Check for local insertion.
                        scalar_t d_sq_loc = tree_local.query(pt_new).second[0];
                        if (d_sq_loc >= (shape.zeta * check_radius) *
                                        (shape.zeta * check_radius)) {
                            local_param_d.addInternalNode(param_new, 1);
                            tree_local.insert(pt_new);
                            end_node++;

                            // Check for global insertion.
                            scalar_t d_sq_global = tree_global.query(pt_new).second[0];
                            if (d_sq_global >= (shape.zeta * check_radius) *
                                               (shape.zeta * check_radius)) {
                                vec_t normal_new = surface_fill_internal::
                                compute_normal(patch.jacobian(param_new, pt_new, w_new));
                                domain.addBoundaryNode(pt_new, type, normal_new);
                                tree_global.insert(pt_new);
                            }
                        }
                    }
                    cur_node++;
                }
            }

            // Random another parameter seed around the middle of the domain, so that fill is
            // successful even if patches overlap.
            std::uniform_real_distribution<> dis(0.0, 1.0);
            param_vec_t param_another = bs.beg() + 0.25 * (bs.end() - bs.beg()) +
                                        dis(gen) * 0.25 * (bs.end() - bs.beg());
            local_param_d.addInternalNode(param_another, 1);

            scalar_t w_another;
            vec_t pt_another = patch.evaluate(param_another, &w_another);
            scalar_t check_radius_another = h(pt_another);
            scalar_t d_sq_another = tree_global.query(pt_another).second[0];
            if (d_sq_another >= (shape.zeta * check_radius_another) *
                                (shape.zeta * check_radius_another)) {
                vec_t normal_another = surface_fill_internal::
                compute_normal(patch.jacobian(param_another, pt_another, w_another));
                domain.addBoundaryNode(pt_another, type, normal_another);
                tree_global.insert(pt_another);
            }

            // Surface fill.
            auto param = [&patch](const param_vec_t& t) {
                return patch.evaluatePointAndJacobian(t);
            };
            gsf.fillParametrization(domain, local_param_d, param, h, tree_global, type);
        }

        return domain;
    }
};

template <typename vec_t>
struct nurbs_shape_internal::NURBSShapeHelper<vec_t, Vec<typename vec_t::scalar_t, 1>>{
    typedef Vec<typename vec_t::scalar_t, 1> param_vec_t;
    typedef typename vec_t::scalar_t scalar_t;
    enum { dim = vec_t::dim,
        proj_dim = dim + 1,
        param_dim = param_vec_t::dim};
    typedef Vec<scalar_t, proj_dim> proj_vec_t;

    static DomainDiscretization<vec_t> discretizeBoundaryWithDensity(
            const NURBSShape<vec_t, param_vec_t>& shape,
            const std::function<scalar_t(vec_t)>& h, int type) {
        if (type == 0) type = -1;
        std::mt19937 gen(shape.seed_);

        GeneralSurfaceFill<vec_t, param_vec_t> gsf;
        gsf.seed(shape.seed_).maxPoints(shape.max_points).proximityTolerance(shape.zeta)
           .numSamples(shape.n_samples);

        // Global domain and KDTree.
        DomainDiscretization<vec_t> domain(shape);
        KDTreeMutable<vec_t> tree_global;

        // Discretize every patch.
        for (const NURBSPatch<vec_t, param_vec_t>& patch : shape.patches) {
            BoxShape<param_vec_t> bs = patch.getDomain();
            DomainDiscretization<param_vec_t> local_param_d(bs);

            // Insert boundary edges and random parameter as seeds.
            std::array<param_vec_t, 2> param_seeds{
                patch.getPatchParameterFromBoundaryParameter({}, 0, shape.epsilon),
                patch.getPatchParameterFromBoundaryParameter({}, 1, shape.epsilon)};
            // Another random seed could be generated here to make this algorithm also work on
            // overlapping curves. This would, however, create two "gaps" of size between h
            // and 2h where the fronts meet, instead of just one.

            for (const Vec1d& param_seed : param_seeds) {
                scalar_t w_seed;
                vec_t pt_seed = patch.evaluate(param_seed, &w_seed);
                local_param_d.addInternalNode(param_seed, 1);

                // Insert into global structure.
                scalar_t check_radius_seed = h(pt_seed);
                scalar_t d_sq = tree_global.size() == 0 ? 10 * check_radius_seed * check_radius_seed
                                                        : tree_global.query(pt_seed).second[0];
                if (d_sq >= (shape.zeta * check_radius_seed) * (shape.zeta * check_radius_seed)) {
                    vec_t normal_seed = surface_fill_internal::
                        compute_normal(patch.jacobian(param_seed, pt_seed, w_seed));
                    domain.addBoundaryNode(pt_seed, type, normal_seed);
                    tree_global.insert(pt_seed);
                }
            }

            // Surface fill.
            auto param = [&patch](const param_vec_t& t) {
                return patch.evaluatePointAndJacobian(t);
            };
            gsf.fillParametrization(domain, local_param_d, param, h, tree_global, type);
        }

        return domain;
    }
};
/// @endcond

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_NURBSSHAPE_HPP_
