#include <medusa/bits/spatial_search/KDTreeMutable.hpp>
#include <medusa/bits/types/Vec.hpp>

/**
 * @file
 * Instantiations of KDTreeMutable.
 */

/// @cond
template class mm::KDTreeMutable<mm::Vec1d>;
template std::ostream& mm::operator<<(std::ostream&, const mm::KDTreeMutable<mm::Vec1d>&);
template class mm::KDTreeMutable<mm::Vec2d>;
template std::ostream& mm::operator<<(std::ostream&, const mm::KDTreeMutable<mm::Vec2d>&);
template class mm::KDTreeMutable<mm::Vec3d>;
template std::ostream& mm::operator<<(std::ostream&, const mm::KDTreeMutable<mm::Vec3d>&);
/// @endcond
