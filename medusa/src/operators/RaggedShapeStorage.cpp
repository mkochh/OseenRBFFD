#include <medusa/bits/operators/RaggedShapeStorage.hpp>
#include <medusa/bits/types/Vec.hpp>
#include <medusa/bits/types/Vec.hpp>

/**
 * @file
 * Instantiations of shape storage for shapes of different sizes.
 */

/// @cond
namespace mm {
template class RaggedShapeStorage<Vec1d>;
template class RaggedShapeStorage<Vec2d>;
template class RaggedShapeStorage<Vec3d>;
}
/// @endcond
