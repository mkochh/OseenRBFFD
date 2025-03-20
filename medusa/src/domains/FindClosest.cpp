#include "medusa/bits/domains/FindClosest.hpp"
#include <medusa/bits/types/Vec.hpp>
#include <medusa/bits/domains/DomainDiscretization.hpp>

/**
 * @file
 * Implementation of FindClosest.
 */

namespace mm {

template void FindClosest::operator()(mm::DomainDiscretization<Vec1d>& domain) const;
template void FindClosest::operator()(mm::DomainDiscretization<Vec2d>& domain) const;
template void FindClosest::operator()(mm::DomainDiscretization<Vec3d>& domain) const;

}  // namespace mm
