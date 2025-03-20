#ifndef MEDUSA_BITS_DOMAINS_FINDBALANCEDSUPPORT_HPP_
#define MEDUSA_BITS_DOMAINS_FINDBALANCEDSUPPORT_HPP_

/**
 * @file
 * Implementation of balanced supports.
 */

#include <medusa/bits/spatial_search/KDTree.hpp>
#include "FindBalancedSupport_fwd.hpp"
#include "Eigen/SVD"

/**
 * @file
 * Implementation of class for finding support nodes.
 */

namespace mm {

/// @cond
template <class vec_t>
Eigen::Matrix<typename vec_t::scalar_t, vec_t::dim, vec_t::dim - 1>
FindBalancedSupport::getFrame(const vec_t& normal)  {
    return normal.transpose().jacobiSvd(Eigen::ComputeFullV).matrixV().rightCols(vec_t::dim-1);
}
/// @endcond

template <typename domain_t>
void mm::FindBalancedSupport::operator()(domain_t& domain) const {
    auto for_which = for_which_;
    if (for_which.empty()) for_which = domain.all();
    auto search_among = search_among_;
    if (search_among.empty()) {
        search_among = Range<int>::seq(domain.size());  // include potential ghost nodes
    }
    assert_msg(!domain.positions().empty(), "Cannot find support in an empty domain.");
    assert_msg(min_support_ > 0, "Support size must be greater than 0, got %d.",
               min_support_);
    assert_msg(min_support_ <= search_among.size(), "Support size (%d) cannot exceed number of "
               "points that we are searching among (%d).", min_support_, search_among.size());
    assert_msg(!for_which.empty(), "Set of nodes for which to find the support is empty.");
    assert_msg(!search_among.empty(), "Set of nodes to search among is empty.");
    for (int x : for_which) {
        assert_msg(0 <= x && x < domain.size(), "Index %d out of range [%d, %d) in for_which.",
                   x, 0, domain.size());
    }
    for (int x : search_among) {
        assert_msg(0 <= x && x < domain.size(), "Index %d out of range [%d, %d) in "
                                                "search_among.", x, 0, domain.size());
    }
    KDTree<typename domain_t::vector_t> tree(domain.positions()[search_among]);
    for (int i : for_which) {
        domain.support(i) = balancedSupport(
                domain, tree, i, min_support_, max_support_, search_among, force_self_);
    }
}

template <class vec_t>
bool FindBalancedSupport::mark_quadrant(
        const vec_t& dx, typename vec_t::scalar_t tol, std::vector<bool>& octants_covered,
        const Eigen::Matrix<typename vec_t::scalar_t, vec_t::dim, -1>& basis) {
    // for each basis vector, determine if it is in the same or opposite direction
    // record that value in `octants`
    int octant = 0;
    Eigen::Matrix<typename vec_t::scalar_t, Eigen::Dynamic, 1> s = basis.transpose()*dx;
    for (int b = 0; b < basis.cols(); ++b) {
        if (std::abs(s[b]) < tol) {
            octant = -1;
            break;
        }
        if (s[b] > 0) octant += (1 << b);
    }
    // mark computed octant as covered
    if (octant != -1 && !octants_covered[octant]) {
        octants_covered[octant] = true;
        return true;
    }
    return false;
}

template <typename domain_t, typename kdtree_t>
Range<int> FindBalancedSupport::balancedSupport(domain_t& domain, const kdtree_t& tree, int i,
                                                int min_support_size, int max_support_size,
                                                const Range<int>& search_among, bool force_self) {
    const int dim = domain_t::dim;
    typedef typename domain_t::vector_t vec_t;
    typedef typename domain_t::scalar_t scalar_t;
    Eigen::Matrix<scalar_t, dim, Eigen::Dynamic> basis =
            Eigen::Matrix<scalar_t, dim, dim>::Identity();
    if (domain.type(i) < 0) basis = getFrame(domain.normal(i));

    int support_size = min_support_size;
    const vec_t& center = domain.pos(i);
    double tol = 1e-6;
    int target = 1 << basis.cols();
    int boundary_neighbours_required = (domain.type(i) < 0) ? 2*(dim-1)+1 : 0;
    std::vector<bool> octants_covered(target, false);
    int pos_idx = 0;
    Range<int> indices;

    while (pos_idx < max_support_size) {
        indices = tree.query(center, support_size).first;
        for (int& j : indices) { j = search_among[j]; }

        if (force_self && indices[0] != i) {
            indices.insert(indices.begin(), i);
        }

        for (; pos_idx < indices.size(); ++pos_idx) {
            vec_t dx = domain.pos(indices[pos_idx]) - center;
            if (domain.type(indices[pos_idx]) < 0 && domain.type(i) < 0) {
                boundary_neighbours_required--;
                target -= mark_quadrant(dx, tol, octants_covered, basis);
            } else if (domain.type(i) > 0) {
                target -= mark_quadrant(dx, tol, octants_covered, basis);
            }
            if (target <= 0 && boundary_neighbours_required <= 0 &&
                pos_idx >= min_support_size-1) {
                return Range<int>(indices.begin(), indices.begin()+pos_idx+1);
            }
        }
        support_size *= 2;
        if (support_size > max_support_size) { support_size = max_support_size; }
    }
    return indices;
}

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_FINDBALANCEDSUPPORT_HPP_
