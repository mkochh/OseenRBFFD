#include "medusa/bits/domains/HalfLinksRefine.hpp"
#include <medusa/bits/types/Vec.hpp>

/**
 * @file
 * Explicit instantiations for common parameters.
 */

// @cond
namespace mm {

template Range<int> HalfLinksRefine::operator()(mm::DomainDiscretization<Vec1d>& domain) const;
template Range<int> HalfLinksRefine::operator()(mm::DomainDiscretization<Vec2d>& domain) const;
template Range<int> HalfLinksRefine::operator()(mm::DomainDiscretization<Vec3d>& domain) const;

}  // namespace mm
// @endcond
