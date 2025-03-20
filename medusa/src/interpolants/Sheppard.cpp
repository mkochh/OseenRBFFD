#include <medusa/bits/interpolants/Sheppard.hpp>
#include <medusa/bits/types/Vec.hpp>

/**
 * @file
 * Instantiation of class for Sheppard's interpolants.
 */

template class mm::SheppardInterpolant<mm::Vec1d, double>;
template class mm::SheppardInterpolant<mm::Vec2d, double>;
template class mm::SheppardInterpolant<mm::Vec3d, double>;
