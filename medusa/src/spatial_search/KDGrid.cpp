#include <medusa/bits/spatial_search/KDGrid.hpp>

/**
 * @file
 * Instantiations of KDGrid.
 */

/// @cond
template class mm::KDGrid<mm::Vec1d>;
template std::ostream& mm::operator<<(std::ostream&, const mm::KDGrid<Vec1d>&);
template class mm::KDGrid<mm::Vec2d>;
template std::ostream& mm::operator<<(std::ostream&, const mm::KDGrid<Vec2d>&);
template class mm::KDGrid<mm::Vec3d>;
template std::ostream& mm::operator<<(std::ostream&, const mm::KDGrid<Vec3d>&);
/// @endcond
