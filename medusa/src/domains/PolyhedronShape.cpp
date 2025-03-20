#include "medusa/bits/domains/PolytopeShape.hpp"
#include "medusa/bits/domains/PolyhedronShape.hpp"
#include "medusa/bits/types/Vec.hpp"

/**
 * @file
 * Instantiation of PolyhedronShape class.
 */

namespace mm {

template class PolyhedronShape<Vec3d>;
template class PolytopeShape<Vec3d>;

}  // namespace mm
