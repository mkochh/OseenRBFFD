#include <medusa/bits/operators/UniformShapeStorage.hpp>
#include <medusa/bits/types/Vec.hpp>

/**
 * @file
 * Instantiations of shape storage for uniform shapes.
 */

/// @cond
namespace mm {
template class UniformShapeStorage<Vec1d>;
template class UniformShapeStorage<Vec2d>;
template class UniformShapeStorage<Vec3d>;
}
/// @endcond
