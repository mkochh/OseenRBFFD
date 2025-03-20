#include "medusa/bits/domains/GeneralSurfaceFill.hpp"

/**
 * @file
 * Instantiations of commonly used GeneralSurfaceFill instances.
 */

template class mm::GeneralSurfaceFill<mm::Vec2d, mm::Vec1d>;
template class mm::GeneralSurfaceFill<mm::Vec3d, mm::Vec2d>;
