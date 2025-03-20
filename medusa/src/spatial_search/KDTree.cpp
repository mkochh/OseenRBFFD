#include <medusa/bits/spatial_search/KDTree.hpp>
#include <medusa/bits/types/Vec.hpp>

/**
 * @file
 * Instantiations of KDTree.
 */

/// @cond
template class mm::KDTree<mm::Vec1d>;
template std::ostream& mm::operator<<(std::ostream&, const mm::KDTree<mm::Vec1d>&);
template class mm::KDTree<mm::Vec2d>;
template std::ostream& mm::operator<<(std::ostream&, const mm::KDTree<mm::Vec2d>&);
template class mm::KDTree<mm::Vec3d>;
template std::ostream& mm::operator<<(std::ostream&, const mm::KDTree<mm::Vec3d>&);
/// @endcond
