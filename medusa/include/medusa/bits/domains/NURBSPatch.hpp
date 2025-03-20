#ifndef MEDUSA_BITS_DOMAINS_NURBSPATCH_HPP_
#define MEDUSA_BITS_DOMAINS_NURBSPATCH_HPP_

/**
 * @file
 * Implementation of the NURBS patch class.
 */

#include "NURBSPatch_fwd.hpp"

#include <medusa/bits/domains/cad_helpers.hpp>
#include <medusa/bits/utils/assert.hpp>

namespace mm {

template <typename vec_t, typename param_vec_t>
NURBSPatch<vec_t, param_vec_t>::NURBSPatch(const NdRange<proj_vec_t>& wcp,
                                           const std::array<Range<scalar_t>, param_dim>& ks,
                                           const std::array<int, param_dim>& in_p) :
                                           weighted_control_points{ wcp }, knots{ ks }, p{ in_p },
                                           der_structure{ nullptr } {
    assert_msg(param_dim <= 2, "Only NURBS curves and surfaces are currently supported.");
}

template <typename vec_t, typename param_vec_t>
NURBSPatch<vec_t, param_vec_t>::NURBSPatch(const NdRange<vec_t>& cp, const NdRange<double>& w,
                                           const std::array<Range<scalar_t>, param_dim>& ks,
                                           const std::array<int, param_dim>& in_p) :
                                           knots{ ks }, p{ in_p }, der_structure{ nullptr } {
    assert_msg(param_dim <= 2, "Only NURBS curves and surfaces are currently supported.");

    weighted_control_points = nurbs_patch_internal::NURBSPatchHelper<vec_t, param_vec_t>::
                                                    cpToWcp(cp, w);
}

template <typename vec_t, typename param_vec_t>
NURBSPatch<vec_t, param_vec_t>::NURBSPatch(const NURBSPatch<vec_t, param_vec_t>& patch) noexcept
        : weighted_control_points{ patch.weighted_control_points },
          knots{ patch.knots }, p { patch.p } {
    if (patch.der_structure != nullptr) {
        der_structure = std::unique_ptr<std::array<NURBSPatch<vec_t, param_vec_t>, param_dim>>(
                new std::array<NURBSPatch<vec_t, param_vec_t>, param_dim>(*patch.der_structure));
    }
}

template <typename vec_t, typename param_vec_t>
NURBSPatch<vec_t, param_vec_t>& NURBSPatch<vec_t, param_vec_t>::operator=(
                                        const NURBSPatch<vec_t, param_vec_t>& other) noexcept {
    weighted_control_points = other.weighted_control_points;
    p = other.p;
    knots = other.knots;

    if (other.der_structure == nullptr && der_structure != nullptr) {
        der_structure = nullptr;
    } else if (other.der_structure != nullptr) {
        der_structure = std::unique_ptr<std::array<NURBSPatch<vec_t, param_vec_t>, param_dim>>(
                new std::array<NURBSPatch<vec_t, param_vec_t>, param_dim>(*other.der_structure));
    }

    return *this;
}

template <typename vec_t, typename param_vec_t>
void NURBSPatch<vec_t, param_vec_t>::computeDerivativeStructure() {
    nurbs_patch_internal::NURBSPatchHelper<vec_t, param_vec_t>::computeDerivativeStructure(*this);
}

template <typename vec_t, typename param_vec_t>
vec_t NURBSPatch<vec_t, param_vec_t>::evaluate(const param_vec_t& t, scalar_t* w) const {
    proj_vec_t weighted_pt = evaluateWeighted(t);
    // Project into the desired dimension.
    vec_t pt;
    for (int i = 0; i < dim; i++) {
        pt(i) = weighted_pt(i) / weighted_pt(proj_dim - 1);
    }
    *w = weighted_pt(proj_dim - 1);
    return pt;
}

template <typename vec_t, typename param_vec_t>
vec_t NURBSPatch<vec_t, param_vec_t>::evaluate(const param_vec_t& t) const {
    scalar_t w;
    return evaluate(t, &w);
}

template <typename vec_t, typename param_vec_t>
typename NURBSPatch<vec_t, param_vec_t>::proj_vec_t
NURBSPatch<vec_t, param_vec_t>::evaluateWeighted(const param_vec_t& t) const {
    return nurbs_patch_internal::NURBSPatchHelper<vec_t, param_vec_t>::evaluateWeighted(*this, t);
}

template <typename vec_t, typename param_vec_t>
Eigen::Matrix<typename NURBSPatch<vec_t, param_vec_t>::scalar_t,
        NURBSPatch<vec_t, param_vec_t>::dim, NURBSPatch<vec_t, param_vec_t>::param_dim>
NURBSPatch<vec_t, param_vec_t>::jacobian(const param_vec_t& t, const vec_t& pt,
                                         const scalar_t& w) const {
    return nurbs_patch_internal::NURBSPatchHelper<vec_t, param_vec_t>::jacobian(*this, t, pt, w);
}

template <typename vec_t, typename param_vec_t>
Eigen::Matrix<typename NURBSPatch<vec_t, param_vec_t>::scalar_t,
        NURBSPatch<vec_t, param_vec_t>::dim, NURBSPatch<vec_t, param_vec_t>::param_dim>
NURBSPatch<vec_t, param_vec_t>::jacobian(const param_vec_t& t) const {
    scalar_t w;
    vec_t pt = evaluate(t, &w);
    return jacobian(t, pt, w);
}

template <typename vec_t, typename param_vec_t>
std::pair<vec_t,
        Eigen::Matrix<
                typename NURBSPatch<vec_t, param_vec_t>::scalar_t,
                NURBSPatch<vec_t, param_vec_t>::dim, NURBSPatch<vec_t, param_vec_t>::param_dim>>
NURBSPatch<vec_t, param_vec_t>::evaluatePointAndJacobian(const param_vec_t& t) const {
    scalar_t w;
    vec_t pt = evaluate(t, &w);
    Eigen::Matrix<scalar_t, dim, param_dim> jm = jacobian(t, pt, w);
    return std::pair<vec_t, Eigen::Matrix<scalar_t, dim, param_dim>>(pt, jm);
}

template <typename vec_t, typename param_vec_t>
BoxShape<param_vec_t> NURBSPatch<vec_t, param_vec_t>::getDomain() const {
    param_vec_t beg, end;
    for (int i = 0; i < param_dim; i++) {
        beg(i) = knots[i].front();
        end(i) = knots[i].back();
    }

    return BoxShape<param_vec_t>(beg, end);
}

template <typename vec_t, typename param_vec_t>
Range<NURBSPatch<vec_t, Vec<typename NURBSPatch<vec_t, param_vec_t>::scalar_t,
    NURBSPatch<vec_t, param_vec_t>::param_dim - 1>>>
    NURBSPatch<vec_t, param_vec_t>::getBoundaries() const {
    return nurbs_patch_internal::NURBSPatchHelper<vec_t, param_vec_t>::getBoundaries(*this);
}

template <typename vec_t, typename param_vec_t>
param_vec_t NURBSPatch<vec_t, param_vec_t>::getPatchParameterFromBoundaryParameter(
        const Vec<scalar_t, param_dim - 1>& t, int i, scalar_t epsilon) const {
    return nurbs_patch_internal::NURBSPatchHelper<vec_t, param_vec_t>::
                                 getPatchParameterFromBoundaryParameter(*this, t, i, epsilon);
}

template <typename vec_t, typename param_vec_t>
Vec<typename NURBSPatch<vec_t, param_vec_t>::scalar_t,
NURBSPatch<vec_t, param_vec_t>::param_dim - 1> NURBSPatch<vec_t, param_vec_t>::
    getBoundaryParameterFromPatchParameter(const param_vec_t& t, int i) const {
    return nurbs_patch_internal::NURBSPatchHelper<vec_t, param_vec_t>::
                                 getBoundaryParameterFromPatchParameter(*this, t, i);
}

/// @cond
template <typename vec_t>
struct nurbs_patch_internal::NURBSPatchHelper<vec_t, Vec<typename vec_t::scalar_t, 2>> {
    typedef Vec<typename vec_t::scalar_t, 2> param_vec_t;
    typedef typename vec_t::scalar_t scalar_t;
    enum {dim = vec_t::dim,
        proj_dim = dim + 1,
        param_dim = param_vec_t::dim};
    typedef Vec<scalar_t, proj_dim> proj_vec_t;

    static Range<Range<proj_vec_t>> cpToWcp(const Range<Range<vec_t>>& cp,
                                     const Range<Range<double>>& w) {
        Range<Range<proj_vec_t>> weighted_control_points;

        weighted_control_points.resize(cp.size());
        for (int i = 0; i < cp.size(); i++) {
            weighted_control_points[i].resize(cp[i].size());
            for (int j = 0; j < cp[i].size(); j++) {
                proj_vec_t temp;

                for (int y = 0; y < dim; y++) {
                    temp(y) = w[i][j] * cp[i][j](y);
                }
                temp(proj_dim - 1) = w[i][j];

                weighted_control_points[i][j] = temp;
            }
        }

        return weighted_control_points;
    }

    static void computeDerivativeStructure(NURBSPatch<vec_t, param_vec_t>& patch) {
        if (patch.der_structure != nullptr) {
            return;
        }

        // Calculate control points for the second dimension.
        Range<Range<proj_vec_t>> v_der_wcp(patch.weighted_control_points.size());
        for (int i = 0; i < patch.weighted_control_points.size(); i++) {
            cad_helpers::generate_b_spline_derivative_control_points(patch.p[1],
                                 patch.weighted_control_points[i], patch.knots[1], v_der_wcp[i]);
        }

        // Calculate knots for the second dimension.
        Range<scalar_t> v_der_knots;
        cad_helpers::generate_b_spline_derivative_knots(patch.knots[1], v_der_knots);

        // Calculate control points for the first dimension.
        // Not efficient, but good enough, since the number of control points is usually
        // sufficiently small. Optimization would require an additional data structure.
        Range<Range<proj_vec_t>> u_der_wcp_initial(patch.weighted_control_points[0].size());
        for (int i = 0; i < patch.weighted_control_points[0].size(); i++) {
            Range<proj_vec_t> temp_range(patch.weighted_control_points.size());
            for (int j = 0; j < patch.weighted_control_points.size(); j++) {
                temp_range[j] = patch.weighted_control_points[j][i];
            }
            cad_helpers::generate_b_spline_derivative_control_points(patch.p[0], temp_range,
                                                             patch.knots[0], u_der_wcp_initial[i]);
        }

        Range<Range<proj_vec_t>> u_der_wcp(patch.weighted_control_points.size() - 1);
        for (int i = 0; i < patch.weighted_control_points.size() - 1; i++) {
            u_der_wcp[i].resize(patch.weighted_control_points[0].size());
            for (int j = 0; j < patch.weighted_control_points[0].size(); j++) {
                u_der_wcp[i][j] = u_der_wcp_initial[j][i];
            }
        }

        // Calculate knots for the first dimension.
        Range<scalar_t> u_der_knots;
        cad_helpers::generate_b_spline_derivative_knots(patch.knots[0], u_der_knots);

        // Construct pair of derivative surfaces.
        patch.der_structure =
                std::unique_ptr<std::array<NURBSPatch<vec_t, param_vec_t>, param_dim>>(
                new std::array<NURBSPatch<vec_t, param_vec_t>, param_dim>{
                         NURBSPatch<vec_t, param_vec_t>(u_der_wcp,
                            std::array<Range<scalar_t>, param_dim>{u_der_knots, patch.knots[1]},
                            std::array<int, param_dim>{patch.p[0] - 1, patch.p[1]}),
                        NURBSPatch<vec_t, param_vec_t>(v_der_wcp,
                            std::array<Range<scalar_t>, param_dim>{patch.knots[0], v_der_knots},
                            std::array<int, param_dim>{patch.p[0], patch.p[1] - 1})});
    }

    static proj_vec_t evaluateWeighted(const NURBSPatch<vec_t, param_vec_t>& patch,
                                const param_vec_t& t) {
        // This can maybe be calculated faster. We are especially limited by evaluating in
        // arbitrary points, since NURBS evaluation can be done much faster on a grid.

        // Evaluate along the second dimension.
        Range<proj_vec_t> secondary_control_points(patch.weighted_control_points.size());
        for (int i = 0; i < patch.weighted_control_points.size(); i++) {
            secondary_control_points[i] = cad_helpers::evaluate_b_spline(t(1), patch.p[1],
                                             patch.weighted_control_points[i], patch.knots[1]);
        }

        // Evaluate along the first dimension.
        return cad_helpers::evaluate_b_spline(t(0), patch.p[0], secondary_control_points,
                                              patch.knots[0]);
    }

    static Eigen::Matrix<scalar_t, dim, param_dim> jacobian(
            const NURBSPatch<vec_t, param_vec_t>& patch, const param_vec_t& t, const vec_t& pt,
            const scalar_t& w) {
        assert_msg(patch.der_structure != nullptr, "Derivative structure not initialised, "
                                           "call initDerivative() before calling this function.");

        // Evaluate point on the first derivative NURBS surface.
        proj_vec_t der_u_pt = (*patch.der_structure)[0].evaluateWeighted(t);

        // Calculate first partial derivative.
        vec_t der_u;
        for (int i = 0; i < dim; i++) {
            der_u(i) = (der_u_pt(i) - der_u_pt(proj_dim - 1) * pt(i)) / w;
        }

        // Evaluate point on the first derivative NURBS surface.
        proj_vec_t der_v_pt = (*patch.der_structure)[1].evaluateWeighted(t);

        // Calculate second partial derivative.
        vec_t der_v;
        for (int i = 0; i < dim; i++) {
            der_v(i) = (der_v_pt(i) - der_v_pt(proj_dim - 1) * pt(i)) / w;
        }

        // Construct Jacobian matrix.
        Eigen::Matrix<scalar_t, dim, 2> jm;

        jm.col(0) << der_u;
        jm.col(1) << der_v;

        return jm;
    }

    static Range<NURBSPatch<vec_t, Vec<scalar_t, param_dim - 1>>> getBoundaries(
            const NURBSPatch<vec_t, param_vec_t>& patch) {
        // Boundaries along the v direction.
        NURBSPatch<vec_t, Vec<scalar_t, param_dim - 1>> b1(patch.weighted_control_points.front(),
                                                           {patch.knots[1]}, {patch.p[1]});
        NURBSPatch<vec_t, Vec<scalar_t, param_dim - 1>> b2(patch.weighted_control_points.back(),
                                                           {patch.knots[1]}, {patch.p[1]});

        // Recalculate the range of control to use along the first dimension.
        Range<proj_vec_t> wcp_front(patch.weighted_control_points.size());
        Range<proj_vec_t> wcp_back(patch.weighted_control_points.size());
        for (int i = 0; i < patch.weighted_control_points.size(); i++) {
            wcp_front[i] = patch.weighted_control_points[i].front();
            wcp_back[i] = patch.weighted_control_points[i].back();
        }

        // Boundaries along the u direction.
        NURBSPatch<vec_t, Vec<scalar_t, param_dim - 1>> b3(wcp_front, {patch.knots[0]},
                                                           {patch.p[0]});
        NURBSPatch<vec_t, Vec<scalar_t, param_dim - 1>> b4(wcp_back, {patch.knots[0]},
                                                           {patch.p[0]});

        return Range<NURBSPatch<vec_t, Vec<scalar_t, param_dim - 1>>>{b1, b2, b3, b4};
    }

    static param_vec_t getPatchParameterFromBoundaryParameter(
            const NURBSPatch<vec_t, param_vec_t>& patch, const Vec<scalar_t, param_dim - 1>& t,
            int i, scalar_t epsilon) {
        scalar_t eps_1 = epsilon * (patch.knots[0].back() - patch.knots[0].front());
        scalar_t eps_2 = epsilon * (patch.knots[1].back() - patch.knots[1].front());

        if (i == 0) return Vec<scalar_t, 2>(patch.knots[0].front() + eps_1, t(0));
        else if (i == 1) return Vec<scalar_t, 2>(patch.knots[0].back() - eps_1, t(0));
        else if (i == 2) return Vec<scalar_t, 2>(t(0), patch.knots[1].front() + eps_2);
        else if (i == 3) return Vec<scalar_t, 2>(t(0), patch.knots[1].back() - eps_2);

        assert_msg(false, "Invalid boundary index, should be 0 <= i <= 3, got %d.", i);
        return {};
    }

    static Vec<scalar_t, param_dim - 1> getBoundaryParameterFromPatchParameter(
            const NURBSPatch<vec_t, param_vec_t>& patch, const param_vec_t& t, int i) {
        if (i == 0 || i == 1) return Vec<scalar_t, param_dim - 1>(t(1));
        else if (i == 2 || i == 3) return Vec<scalar_t, param_dim - 1>(t(0));

        assert_msg(false, "Invalid edge index, should be 0 <= i <= 3, got %d.", i);
        return {};
    }
};

template <typename vec_t>
struct nurbs_patch_internal::NURBSPatchHelper<vec_t, Vec<typename vec_t::scalar_t, 1>> {
    typedef Vec<typename vec_t::scalar_t, 1> param_vec_t;
    typedef typename vec_t::scalar_t scalar_t;
    enum {dim = vec_t::dim,
        proj_dim = dim + 1,
        param_dim = param_vec_t::dim};
    typedef Vec<scalar_t, proj_dim> proj_vec_t;

    static Range<proj_vec_t> cpToWcp(const Range<vec_t>& cp, const Range<double>& w) {
        Range<proj_vec_t> weighted_control_points;

        weighted_control_points.resize(cp.size());
        for (int i = 0; i < cp.size(); i++) {
            proj_vec_t temp;

            for (int j = 0; j < dim; j++) {
                temp(j) = w[i] * cp[i](j);
            }
            temp(proj_dim - 1) = w[i];

            weighted_control_points[i] = temp;
        }

        return weighted_control_points;
    }

    static void computeDerivativeStructure(NURBSPatch<vec_t, param_vec_t>& patch) {
        if (patch.der_structure != nullptr || patch.p[0] == 0) {
            return;
        }

        Range<proj_vec_t> derivative_weighted_control_points;
        Range<scalar_t> derivative_knots;
        cad_helpers::generate_b_spline_derivative(patch.p[0], patch.weighted_control_points,
                          patch.knots[0], derivative_weighted_control_points, derivative_knots);

        patch.der_structure =
                std::unique_ptr<std::array<NURBSPatch<vec_t, param_vec_t>, param_dim>>(
                new std::array<NURBSPatch<vec_t, param_vec_t>, param_dim>{
                    NURBSPatch<vec_t, param_vec_t>(derivative_weighted_control_points,
                                                   {derivative_knots}, {patch.p[0] - 1})});
    }

    static proj_vec_t evaluateWeighted(const NURBSPatch<vec_t, param_vec_t>& patch,
                                const param_vec_t& t) {
        return cad_helpers::evaluate_b_spline(t(0), patch.p[0], patch.weighted_control_points,
                                              patch.knots[0]);
    }

    static Eigen::Matrix<scalar_t, dim, param_dim> jacobian(
            const NURBSPatch<vec_t, param_vec_t>& patch, const param_vec_t& t, const vec_t& pt,
            const scalar_t& w) {
        if (patch.p[0] == 0) {
            return vec_t();
        }
        assert_msg(patch.der_structure != nullptr, "Derivative structure not initialised, "
                                         "call initDerivative() before calling this function.");

        // Evaluate point on the derivative NURBS.
        proj_vec_t der_pt = (*patch.der_structure)[0].evaluateWeighted(t);

        // Calculate Jacobian
        vec_t jacobian;
        for (int i = 0; i < dim; i++) {
            jacobian(i) = (der_pt(i) - der_pt(proj_dim - 1) * pt(i)) / w;
        }

        return jacobian;
    }

    static Range<NURBSPatch<vec_t, Vec<scalar_t, param_dim - 1>>> getBoundaries(
            const NURBSPatch<vec_t, param_vec_t>& patch) {
        return {NURBSPatch<vec_t, Vec<scalar_t, param_dim - 1>>({}, {}, {}),
                NURBSPatch<vec_t, Vec<scalar_t, param_dim - 1>>({}, {}, {})};
    }

    static param_vec_t getPatchParameterFromBoundaryParameter(
            const NURBSPatch<vec_t, param_vec_t>& patch, const Vec<scalar_t, param_dim - 1>& t,
            int i, scalar_t epsilon) {
        BoxShape<param_vec_t> dom = patch.getDomain();
        scalar_t eps = epsilon * (patch.knots[0].back() - patch.knots[0].front());

        if (i == 0) return dom.beg() + param_vec_t(eps);
        else if (i == 1) return dom.end() - param_vec_t(eps);

        assert_msg(false, "Invalid boundary index, should be 0 <= i <= 1, got %d.", i);
        return {};
    }

    static Vec<scalar_t, param_dim - 1> getBoundaryParameterFromPatchParameter(
            const NURBSPatch<vec_t, param_vec_t>& patch, const param_vec_t& t, int i) {
        if (i == 0 || i == 1) return {};

        assert_msg(false, "Invalid edge index, should be 0 <= i <= 1, got %d.", i);
        return {};
    }
};
/// @endcond

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_NURBSPATCH_HPP_
