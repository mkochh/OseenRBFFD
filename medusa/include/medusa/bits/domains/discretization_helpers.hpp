#ifndef MEDUSA_BITS_DOMAINS_DISCRETIZATION_HELPERS_HPP_
#define MEDUSA_BITS_DOMAINS_DISCRETIZATION_HELPERS_HPP_

/**
 * @file
 * Contains helpers for discretizing various geometric objects.
 *
 * @example test/domains/discretization_helpers_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/types/Range.hpp>
#include <Eigen/Core>
#include <random>

namespace mm {

/**
 * Namespace containing helpers for discretizing various geometric objects.
 */
namespace discretization_helpers {

/**
 * Returns nodes lying on the line segment `pq` with approximate distances `delta_r`.
 * Boundary point are not included in return set, but are included in the spacing.
 * The last discretization point almost never falls on point q, but before it, and
 * we decide on whether to include it or not, based on whether spacing at `q` is more
 * violated with including the point or without.
 * @return Discretization points spaced according to the `delta_r` with last and first point
 * (i.e.\ `p` and `q`) excluded.
 * @ingroup domains
 * @snippet domains/discretization_helpers_test.cpp discretizeLineWithDensity
 */
template <typename vec_t, typename func_t>
Range<vec_t> discretizeLineWithDensity(const vec_t& p, const vec_t& q, const func_t& delta_r) {
    typedef typename vec_t::scalar_t scalar_t;
    scalar_t dist = (q-p).norm();
    assert_msg(dist >= 1e-15, "Points %s and %s are equal.", p, q);
    vec_t dir = (q-p) / dist;
    scalar_t cur_dist = 0;
    Range<scalar_t> dists;
    while (cur_dist < dist) {
        cur_dist += delta_r(p+cur_dist*dir);
        dists.push_back(cur_dist);
    }
    dists.pop_back();
    if (dists.size() <= 2) return {};
    // decide whether to include the last point
    scalar_t spacing_with_last_point = dist - dists.back();
    scalar_t spacing_without_last_point = dist - dists[dists.size()-2];
    scalar_t desired_spacing = delta_r(q);
    scalar_t error_with = std::abs(1 - spacing_with_last_point / desired_spacing);
    scalar_t error_without = std::abs(1 - spacing_without_last_point / desired_spacing);
    if (error_without < error_with) {
        dists.pop_back();
    }

    Range<vec_t> points(dists.size());
    for (int i = 0; i < dists.size(); ++i) {
        points[i] = p + dists[i]*dir;
    }
    return points;
}

/**
 * Discretizes a sphere with given radius uniformly with `num_points` points on the great circle.
 * This is a class because C++ does not allow partial template specializations of functions.
 * @tparam dim Dimension of the sphere.
 * @tparam scalar_t Data type for numeric computations, e.g.\ `double` or `float`.
 *
 * @snippet domains/discretization_helpers_test.cpp SphereDiscretization usage example
 * @ingroup domains
 */
template <typename scalar_t, int dim>
struct SphereDiscretization {
    /**
     * Construct a randomized discretization.
     * @param radius Radius of the sphere.
     * @param num_samples Number of points on the equator, implies nodal spacing `dp = 2*pi*r/n`.
     * @param generator A random number generator.
     * @return A vector of discretization points.
     */
    template <typename generator_t>
    static std::vector<Eigen::Matrix<scalar_t, dim, 1>,
                       Eigen::aligned_allocator<Eigen::Matrix<scalar_t, dim, 1>>> construct(
            scalar_t radius, int num_samples, generator_t& generator) {
        scalar_t dphi = 2 * Pi<scalar_t>::value / num_samples;
        std::uniform_real_distribution<scalar_t> distribution(0, Pi<scalar_t>::value);
        scalar_t offset = distribution(generator);
        std::vector<Eigen::Matrix<scalar_t, dim, 1>,
                Eigen::aligned_allocator<Eigen::Matrix<scalar_t, dim, 1>>> result;
        for (int i = 0; i < num_samples / 2; ++i) {
            scalar_t phi = i * dphi + offset;
            if (phi > Pi<scalar_t>::value) phi -= Pi<scalar_t>::value;
            int slice_n = static_cast<int>(std::ceil(num_samples * std::sin(phi)));
            if (slice_n == 0) continue;
            auto slice = SphereDiscretization<scalar_t, dim - 1>::construct(
                    radius * std::sin(phi), slice_n, generator);
            Eigen::Matrix<scalar_t, dim, 1> v;
            for (const auto& p : slice) {
                v[0] = radius * std::cos(phi);
                v.template tail<dim - 1>() = p;
                result.push_back(v);
            }
        }
        return result;
    }
    /// Construct the discretization.
    static std::vector<Eigen::Matrix<scalar_t, dim, 1>,
            Eigen::aligned_allocator<Eigen::Matrix<scalar_t, dim, 1>>> construct(
            scalar_t radius, int num_samples) {
        scalar_t dphi = 2 * Pi<scalar_t>::value / num_samples;
        std::vector<Eigen::Matrix<scalar_t, dim, 1>,
                Eigen::aligned_allocator<Eigen::Matrix<scalar_t, dim, 1>>> result;
        for (int i = 0; i < num_samples / 2; ++i) {
            scalar_t phi = i * dphi;
            if (phi > Pi<scalar_t>::value) phi -= Pi<scalar_t>::value;
            int slice_n = static_cast<int>(std::ceil(num_samples * std::sin(phi)));
            if (slice_n == 0) continue;
            auto slice = SphereDiscretization<scalar_t, dim - 1>::construct(
                    radius * std::sin(phi), slice_n);
            Eigen::Matrix<scalar_t, dim, 1> v;
            for (const auto& p : slice) {
                v[0] = radius * std::cos(phi);
                v.template tail<dim - 1>() = p;
                result.push_back(v);
            }
        }
        return result;
    }
};

/// Two-dimensional base case of the discretization.
template <typename scalar_t>
struct SphereDiscretization<scalar_t, 2> {
    /// Construct a randomized discretization.
    template <typename generator_t>
    static std::vector<Eigen::Matrix<scalar_t, 2, 1>,
            Eigen::aligned_allocator<Eigen::Matrix<scalar_t, 2, 1>>> construct(
            scalar_t radius, int num_samples, generator_t& generator) {
        scalar_t dphi = 2 * Pi<scalar_t>::value / num_samples;
        std::uniform_real_distribution<scalar_t> distribution(0, 2 * Pi<scalar_t>::value);
        scalar_t offset = distribution(generator);
        std::vector<Eigen::Matrix<scalar_t, 2, 1>,
                Eigen::aligned_allocator<Eigen::Matrix<scalar_t, 2, 1>>> result;
        for (int i = 0; i < num_samples; ++i) {
            scalar_t phi = i * dphi + offset;
            result.emplace_back(radius * std::cos(phi), radius * std::sin(phi));
        }
        return result;
    }
    /// Construct the discretization.
    static std::vector<Eigen::Matrix<scalar_t, 2, 1>,
            Eigen::aligned_allocator<Eigen::Matrix<scalar_t, 2, 1>>> construct(
            scalar_t radius, int num_samples) {
        scalar_t dphi = 2 * Pi<scalar_t>::value / num_samples;
        std::vector<Eigen::Matrix<scalar_t, 2, 1>,
                Eigen::aligned_allocator<Eigen::Matrix<scalar_t, 2, 1>>> result;
        for (int i = 0; i < num_samples; ++i) {
            scalar_t phi = i * dphi;
            result.emplace_back(radius * std::cos(phi), radius * std::sin(phi));
        }
        return result;
    }
};

/// One-dimensional base case of the discretization.
template <typename scalar_t>
struct SphereDiscretization<scalar_t, 1> {
    template <typename generator_t>
    /// Construct a randomized discretization.
    static std::vector<Eigen::Matrix<scalar_t, 1, 1>,
            Eigen::aligned_allocator<Eigen::Matrix<scalar_t, 1, 1>>> construct(
            scalar_t radius, int, generator_t&) {
        return {Eigen::Matrix<scalar_t, 1, 1>(-radius), Eigen::Matrix<scalar_t, 1, 1>(radius)};
    }
    /// Construct the discretization.
    static std::vector<Eigen::Matrix<scalar_t, 1, 1>,
            Eigen::aligned_allocator<Eigen::Matrix<scalar_t, 1, 1>>> construct(
            scalar_t radius, int) {
        return {Eigen::Matrix<scalar_t, 1, 1>(-radius), Eigen::Matrix<scalar_t, 1, 1>(radius)};
    }
};

}  // namespace discretization_helpers
}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_DISCRETIZATION_HELPERS_HPP_
