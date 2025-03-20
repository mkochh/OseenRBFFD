#ifndef MEDUSA_BITS_DOMAINS_GRAINDROPFILL_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_GRAINDROPFILL_FWD_HPP_

/**
 * @file
 * Declaration of the grain drop node placing algorithm.
 *
 * @example test/domains/GrainDropFill_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/types/Vec_fwd.hpp>

/**
 * @file
 * Declaration of GrainDropFill class.
 */

namespace mm {

template <typename vec_t>
class DomainDiscretization;

/**
 * Implements general `n`-d node placing algorithm, as described in https://arxiv.org/abs/1906.00636
 * This node generation algorithm simulates dropping grains with variable radius in a box.
 * First, a large bounding box is filled and only the nodes in the domain shape are kept.
 * Thus, gaps may appear at the boundary and BasicRelax can be used to smooth them.
 *
 * Usage example:
 * @snippet domains/GrainDropFill_test.cpp GDF usage example
 *
 * @sa GeneralFill, BasicRelax
 * @ingroup domains
 */
template <typename vec_t>
class GrainDropFill {
  public:
    typedef typename vec_t::scalar_t scalar_t;  ///< Scalar type;
    typedef vec_t vector_t;  ///< Vector type.
    /// Store dimension of the domain.
    enum { /** Dimensionality of the domain. */ dim = vec_t::dim };

  private:
    int max_points = 5000000;  ///< Maximal number of points generated.
    vec_t bot;  ///< Bottom left corner of the box to fill.
    vec_t top;  ///< Top right corner of the box to fill.
    scalar_t dx = -1.0;  ///< Initial discretization spacing.
    int seed_;  ///<  Seed for the random number generator.
    scalar_t excess_factor = 10;  ///< Fill `excess_factor*dx` over the top of the domain.

  public:
    /// Prepare to fill domain within `[bot, top]` bounding box.
    GrainDropFill(const vec_t& bot, const vec_t& top);

    /// Maximal number of points generated.
    GrainDropFill& maxPoints(int max_points) { this->max_points = max_points; return *this; }
    /// Set custom seed for the random number generator.
    GrainDropFill& seed(int seed) { seed_ = seed; return *this; }
    /**
     * Set initial discretization density by specifying number of dots along each dimension.
     * The recommended density is 10 times less than the minimal nodal spacing.
     * Number of data points along the last dimension is ignored.
     */
    GrainDropFill& initialSpacing(scalar_t dx) { this->dx = dx; return *this; }
    /// Set percentage of allowed generation over the top.
    GrainDropFill& excessFactor(scalar_t factor) { excess_factor = factor; return *this; }

    /**
     * Fills domain with a quality node distribution.
     * @param domain Domain to fill with nodes. Can be partially filled.
     * @param h Nodal spacing function.
     * @param type Type of the nodes. If 0 the engines default value is used.
     *
     * @warning Header <code> \#include <medusa/bits/domains/GrainDropFill.hpp></code> must be
     * included additionally to <code> \#include <medusa/Medusa_fwd.hpp></code>.
     */
    template <typename func_t>
    void operator()(DomainDiscretization<vec_t>& domain, const func_t& h, int type = 0) const;

    /// Overload for constant function.
    void operator()(DomainDiscretization<vec_t>& domain, const scalar_t& h, int type = 0) const {
        this->operator()(domain, [=] (const vector_t&) { return h; }, type);
    }

    /// Actual fill algorithm which fills `[bot, top]` box with spacing `h`.
    template <typename func_t>
    std::vector<vec_t> fillBox(const func_t& h) const;
};

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_GRAINDROPFILL_FWD_HPP_
