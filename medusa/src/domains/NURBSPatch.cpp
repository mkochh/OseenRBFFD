#include "medusa/bits/domains/NURBSPatch.hpp"

/**
 * @file
 * Instantiations of commonly used NURBSPatch instances.
 */

/// @cond
template class mm::NURBSPatch<mm::Vec3d, mm::Vec2d>;

template class mm::NURBSPatch<mm::Vec3d, mm::Vec1d>;
template class mm::NURBSPatch<mm::Vec2d, mm::Vec1d>;
/// @endcond
