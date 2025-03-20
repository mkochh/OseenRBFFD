#ifndef MEDUSA_BITS_DOMAINS_POLYTOPESHAPE_FWD_HPP_
#define MEDUSA_BITS_DOMAINS_POLYTOPESHAPE_FWD_HPP_

#include <medusa/Config.hpp>
#include "DomainShape_fwd.hpp"
#include "PolygonShape_fwd.hpp"

/**
 * @file
 * Declaration of PolytopeShape class.
 *
 * @example test/domains/PolytopeShape_test.cpp
 */

namespace mm {

/**
 * Shape representing a simple nonempty polytope (i.e., non-self intersecting polygon in
 * 2D and polyhedron in 3D.). Can be used as a generic dimension independent
 * substitute for PolygonShape or PolyhedronShape.
 *
 * @tparam vec_t This shape class is used only for 2D and 3D domains. General polytopes are not
 * supported.
 *
 * @warning 3D instantiation of this class depends on CGAL and is not included by default. You must
 * explicitly include PolyhedronShape.hpp.
 *
 * Usage example:
 * @snippet domains/PolytopeShape_test.cpp PolytopeShape usage example
 *
 * @ingroup domains
 * @sa PolygonShape, PolyhedronShape
 */
template <typename vec_t>
class PolytopeShape : public DomainShape<vec_t> {
    static_assert(vec_t::dim == 2 || vec_t::dim == 3, "Only available in 2 and 3 dimensions");
};

// This must be hidden during a doxygen run due to a bug in processing const declarations with
// differently named template parameters.
// See https://github.com/doxygen/doxygen/issues/8178
#ifndef DOXYGEN_RUN
/**
 * Specialization for 2D, implementing polygons.
 * @sa PolygonShape
 * @ingroup domains
 */
template <typename Scalar>
class PolytopeShape<Vec<Scalar, 2>> : public PolygonShape<Vec<Scalar, 2>> {
    using base_t = PolygonShape<Vec<Scalar, 2>>;
    using typename base_t::PolygonShape;
  public:
    PolytopeShape(const base_t& shape) : base_t(shape) {}
};
#endif

// This must be hidden during a doxygen run due to a bug in processing forward declarations.
// See https://github.com/doxygen/doxygen/issues/8177
#ifndef DOXYGEN_RUN
template <typename Scalar>
class PolytopeShape<Vec<Scalar, 3>>;   // You must include PolyhedronShape.hpp for the definition.
#endif

}  // namespace mm

#endif  // MEDUSA_BITS_DOMAINS_POLYTOPESHAPE_FWD_HPP_
