#include <medusa/bits/domains/DomainDiscretization.hpp>
#include <medusa/bits/types/Vec.hpp>

/**
 * @file
 * Instantiations for common domain discretizations.
 */

template class mm::DomainDiscretization<mm::Vec1d>;
template mm::Range<int> mm::DomainDiscretization<mm::Vec1d>::reorderNodes(const std::less<std::pair<mm::Vec1d, int>>& cmp);  // NOLINT
template class mm::DomainDiscretization<mm::Vec2d>;
template mm::Range<int> mm::DomainDiscretization<mm::Vec2d>::reorderNodes(const std::less<std::pair<mm::Vec2d, int>>& cmp);  // NOLINT
template class mm::DomainDiscretization<mm::Vec3d>;
template mm::Range<int> mm::DomainDiscretization<mm::Vec3d>::reorderNodes(const std::less<std::pair<mm::Vec3d, int>>& cmp);  // NOLINT
