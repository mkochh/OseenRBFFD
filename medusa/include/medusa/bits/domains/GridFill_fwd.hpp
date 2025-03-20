#ifndef MEDUSA_BITS_DOMAINS_GRIDFILL_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_GRIDFILL_FWD_HPP_

/**
 * @file
 * Declaration of the GridFill node generation algorithm.
 *
 * @example test/domains/GridFill_test.cpp
 */

#include <medusa/Config.hpp>

namespace mm {

template <typename vec_t>
class DomainDiscretization;

/**
 * Discretizes the domain using `n`-dimensional grid with given endpoints and spacing `h`.
 * Only nodes that are inside the domain and at least `h` away from existing nodes are added.
 * Gaps may appear at the boundary and BasicRelax can be used to smooth them.
 *
 * Usage example:
 * @snippet domains/GridFill_test.cpp GridFill usage example
 *
 * @sa GeneralFill, BasicRelax
 * @ingroup domains
 */
template <typename vec_t>
class GridFill {
  public:
    typedef typename vec_t::scalar_t scalar_t;  ///< Scalar type;
    typedef vec_t vector_t;  ///< Vector type.
    /// Store dimension of the domain.
    enum { /** Dimensionality of the domain. */ dim = vec_t::dim };

  private:
    vec_t bot;  ///< Bottom left corner of the box to fill.
    vec_t top;  ///< Top right corner of the box to fill.

  public:
    /// Prepare to fill domain within `[bot, top]` bounding box.
    GridFill(const vec_t& bot, const vec_t& top);

    /**
     * Fills domain with a grid node distribution with spacing `h`.
     * Newly added nodes are at least `h` away from all previous nodes.
     * @param domain Domain to fill with nodes. Can be partially filled.
     * @param h Nodal spacing.
     * @param type Type of the nodes. If 0 the engines default value is used.
     * @warning The resulting distribution can be rough near the boundary. Consider using relaxation
     * techniques or different fill engines, if this is a problem.
     */
    void operator()(DomainDiscretization<vec_t>& domain, const scalar_t& h, int type = 0) const;
};


}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_GRIDFILL_FWD_HPP_
