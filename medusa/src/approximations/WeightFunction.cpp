#include <medusa/bits/approximations/WeightFunction.hpp>
#include <medusa/bits/approximations/Gaussian.hpp>

/**
 * @file
 * Explicit instantiations for common parameters.
 */

#include <medusa/bits/types/Vec.hpp>

template class mm::RBFWeight<mm::Gaussian<double>, mm::Vec1d>;
template class mm::RBFWeight<mm::Gaussian<double>, mm::Vec2d>;
template class mm::RBFWeight<mm::Gaussian<double>, mm::Vec3d>;
template class mm::NoWeight<mm::Vec1d>;
template class mm::NoWeight<mm::Vec2d>;
template class mm::NoWeight<mm::Vec3d>;
