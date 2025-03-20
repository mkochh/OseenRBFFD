#ifndef MEDUSA_BITS_DOMAINS_GRAINDROPFILL_HPP_
#define MEDUSA_BITS_DOMAINS_GRAINDROPFILL_HPP_

/**
 * @file
 * Implementation of the GrainDropFill fill engine.
 */

#include "GrainDropFill_fwd.hpp"
#include <medusa/bits/utils/randutils.hpp>
#include <medusa/bits/spatial_search/KDTree.hpp>
#include <medusa/bits/spatial_search/Grid.hpp>
#include <medusa/bits/domains/discretization_helpers.hpp>

namespace mm {

template <typename vec_t>
GrainDropFill<vec_t>::GrainDropFill(const vec_t& bot, const vec_t& top) :
        bot(bot), top(top), seed_(get_seed()) {
    for (int i = 0; i < dim; ++i) {
        assert_msg(top[i] > bot[i], "Bottom must be less than top, got bot = %s, top = %s.",
                   bot, top);
    }
}

template <typename vec_t>
template <typename func_t>
void GrainDropFill<vec_t>::operator()(DomainDiscretization<vec_t>& domain, const func_t& h,
                                      int type) const {
    if (type == 0) type = 1;
    std::vector<vec_t> nodes = fillBox(h);
    KDTree<vec_t> boundary_tree(domain.positions());
    for (const vec_t& p : nodes) {
        if (!domain.contains(p, boundary_tree)) continue;
        scalar_t r2 = boundary_tree.query(p).second[0];
        scalar_t allowed_r = 0.75*h(p);
        if (r2 >= allowed_r*allowed_r) {
            domain.addInternalNode(p, type);
        }
    }
}

template <typename vec_t>
template <typename func_t>
std::vector<vec_t> GrainDropFill<vec_t>::fillBox(const func_t& h) const {
    assert_msg(dx > 0, "Initial spacing must be positive, got %g. Did you set it correctly?", dx);

    std::vector<vec_t> result;

    // Initialization.
    Vec<int, dim-1> grid_size;
    for (int i = 0; i < dim-1; ++i) {
        grid_size[i] = iceil((top[i] - bot[i]) / dx) + 1;
        assert_msg(grid_size[i] >= 1, "Grid size in dimension %d is not positive.");
    }

    Grid<scalar_t, dim-1, int, Vec<int, dim-1>> z_coor(grid_size, bot[dim-1]);
    Eigen::Map<Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>>
    eig_z(z_coor.data().data(), z_coor.size());
    eig_z += 1e-4*dx*Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>::Random(z_coor.size());

    // Main iteration loop. Find locally minimal point and raise its neighbourhood up.
    int min_lin_idx;
    scalar_t min_z = eig_z.minCoeff(&min_lin_idx);
    Vec<int, dim-1> min_idx = z_coor.multiIndex(min_lin_idx);
    vec_t p;
    int dot_nr = 0;
    while (dot_nr < max_points && min_z <= top[dim-1] + excess_factor*dx) {
        // Create new point.
        for (int i = 0; i < dim-1; ++i) p[i] = bot[i] + min_idx[i] * dx;
        p[dim-1] = min_z;
        result.push_back(p);
        ++dot_nr;

        // Update heights of neighbours.
        scalar_t r = h(p);
        int idx_rad = ifloor(r / dx);
        Vec<int, dim-1> low = (min_idx - Vec<int, dim-1>::Constant(idx_rad)).cwiseMax(0);
        Vec<int, dim-1> high = (min_idx + Vec<int, dim-1>::Constant(idx_rad+1)).cwiseMin(grid_size);
        Vec<int, dim-1> cnt = low;
        int cnt_lin;
        do {
            scalar_t height_update = r*r-((cnt-min_idx).template cast<scalar_t>()*dx).squaredNorm();
            cnt_lin = z_coor.linearIndex(cnt);
            if (height_update > 0) {
                z_coor[cnt_lin] = std::max(z_coor[cnt_lin], min_z + std::sqrt(height_update));
            }
        } while (incrementCounter(cnt, low, high));
        min_z = z_coor[min_lin_idx];  // Update minimum to be able to find something smaller later.

        // Find next local minimum
        int dist = 2*idx_rad;
        while (dist > idx_rad) {
            for (int i = 0; i < dim-1; ++i) {
                low[i] = (((min_idx[i] - 2*idx_rad) % grid_size[i]) + grid_size[i]) % grid_size[i];
                high[i] = (min_idx[i] + 2*idx_rad + 1) % grid_size[i];
            }
            auto old_min_idx = min_idx;
            cnt = low;
            do {
                cnt_lin = z_coor.linearIndex(cnt);
                if (z_coor[cnt_lin] < min_z) {
                    min_z = z_coor[cnt_lin];
                    min_lin_idx = cnt_lin;
                    min_idx = cnt;
                }
            } while (incrementCyclicCounter(cnt, low, high, grid_size));
            dist = (min_idx - old_min_idx).cwiseAbs().maxCoeff();
        }
    }
    return result;
}

/// Specialization for 1D
template <>
template <typename func_t>
std::vector<Vec1d> GrainDropFill<Vec1d>::fillBox(const func_t& h) const {
    auto result = discretization_helpers::discretizeLineWithDensity(bot, top, h);
    result.resize(std::min(result.size(), max_points));
    return result;
}

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_GRAINDROPFILL_HPP_
