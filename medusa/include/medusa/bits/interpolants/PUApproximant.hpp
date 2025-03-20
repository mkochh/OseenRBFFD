
#ifndef MEDUSA_BITS_INTERPOLANTS_PUAPPROXIMANT_HPP_
#define MEDUSA_BITS_INTERPOLANTS_PUAPPROXIMANT_HPP_

/**
 * @file
 * Implementation of partition of unity approximant class.
 */

#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/spatial_search/KDTree.hpp>
#include <medusa/bits/utils/numutils.hpp>

#include "PUApproximant_fwd.hpp"

namespace mm {

/// @cond
template <typename vec_t>
template <typename scalar_t, typename engine_t>
Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> PUApproximant<vec_t>::evaluate(
    const DomainDiscretization<vec_t>& domain, const eigen_vt<scalar_t>& values,
    const Range<vec_t>& query_points, d_scalar_t radius_factor, const engine_t& engine) {
    KDTree<vec_t> tree(query_points);  // Build KDTree.

    return evaluate(domain, values, tree, radius_factor, engine);
}

template <typename vec_t>
template <typename scalar_t, typename engine_t>
Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> PUApproximant<vec_t>::evaluate(
    const DomainDiscretization<vec_t>& domain, const eigen_vt<scalar_t>& values,
    const KDTree<vec_t>& query_points_tree, d_scalar_t radius_factor, const engine_t& engine) {
    int N = domain.size();
    assert_msg(N >= 0, "Missing domain positions.");
    for (int i = 0; i < domain.size(); i++) {
        int support_size = domain.supportSize(i);
        if (support_size <= 0) {
            assert_msg(support_size >= 0,
                       "Expected support size >= 0. Node %d has support size %d. Did you forget "
                       "to call findSupport()?",
                       i, support_size);
        }
    }
    int M = values.size();
    assert_msg(M == N,
               "Number of function values given different from number of "
               "nodes in the domain. Got %d function values and %d nodes.",
               M, N);
    assert_msg(radius_factor >= 0 && radius_factor <= 1,
               "Effective radius must be "
               "larger than 0 and smaller than 1, got %g",
               radius_factor);

    int Nq = query_points_tree.size();  // Number of query points.
    eigen_vt<scalar_t> result_num(Nq);
    result_num.setConstant(0);
    eigen_vt<d_scalar_t> result_den(Nq);
    result_den.setConstant(0);
    Eigen::VectorXi count(Nq);
    count.setConstant(0);

#pragma omp declare reduction(+: eigen_vt<scalar_t>: omp_out = omp_out + omp_in)\
        initializer(omp_priv = eigen_vt<scalar_t>::Zero(omp_orig.size()))
#pragma omp declare reduction(+: Eigen::VectorXi: omp_out = omp_out + omp_in)\
        initializer(omp_priv = Eigen::VectorXi::Zero(omp_orig.size()))

#pragma omp parallel for reduction(+ : result_num, result_den, count)
    for (int i = 0; i < N; ++i) {
        if (domain.supportSize(i) <= 0) continue;  // Check i-th support size.

        const vec_t& p = domain.pos(i);             // Central node.
        Range<vec_t> sup = domain.supportNodes(i);  // Support nodes.
        auto dr = (sup[0] - sup[sup.size() - 1]).norm();
        d_scalar_t radius = radius_factor * dr;
        d_scalar_t radius2 = ipow<2>(radius);  // Query radius.

        Range<int> idx;
        Range<double> dist2;
        std::tie(idx, dist2) =
            query_points_tree.query(p, radius2);  // Query for close enough nodes.
        int n = idx.size();
        if (n == 0) continue;

        // Get approximants.
        auto appr = engine.getApproximant(p, domain.supportNodes(i), values(domain.support(i)));

        for (int j = 0; j < n; ++j) {
            const int& k = idx[j];
            ++count[k];
            d_scalar_t w = weight(std::sqrt(dist2[j] / radius2));
            result_num[k] += w * appr(query_points_tree.get(k));
            result_den[k] += w;
        }
    }

    // Check if all query nodes are covered.
    for (int i = 0; i < Nq; ++i) {
        assert_msg(count[i] > 0,
                   "Query point %d has an undefined value. Point %s is not covered by any of "
                   "the balls with centers in discretization nodes. Consider increasing the "
                   "ball radius.",
                   i, query_points_tree.get(i));
    }

    return result_num.cwiseQuotient(result_den);
}
/// @endcond
}  // namespace mm

#endif  // MEDUSA_BITS_INTERPOLANTS_PUAPPROXIMANT_HPP_
