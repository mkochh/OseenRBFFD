#include <medusa/bits/domains/GeneralFill.hpp>
#include <medusa/bits/domains/DomainDiscretization_fwd.hpp>
#include <medusa/bits/types/Vec.hpp>

/**
 * @file
 * Instantiations of Poisson disk sampling algorithm.
 */

template class mm::GeneralFill<mm::Vec1d>;
template class mm::GeneralFill<mm::Vec2d>;
template class mm::GeneralFill<mm::Vec3d>;
