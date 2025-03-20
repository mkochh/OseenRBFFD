#include "medusa/bits/domains/NURBSShape.hpp"

/**
 * @file
 * Instantiations of commonly used NURBSShape instances.
 */

/// @cond
template class mm::NURBSShape<mm::Vec3d, mm::Vec2d>;

template class mm::NURBSShape<mm::Vec3d, mm::Vec1d>;
template class mm::NURBSShape<mm::Vec2d, mm::Vec1d>;
/// @endcond
