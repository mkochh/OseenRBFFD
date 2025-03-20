#include <medusa/bits/domains/BoxShape.hpp>
#include <medusa/bits/types/Vec.hpp>

/**
 * @file
 * Instantiation of class for box shaped domains.
 */

template class mm::BoxShape<mm::Vec1d>;
template class mm::BoxShape<mm::Vec2d>;
template class mm::BoxShape<mm::Vec3d>;
