#include "medusa/bits/approximations/RBFBasis.hpp"
#include <medusa/bits/types/Vec.hpp>
#include <medusa/bits/approximations/Gaussian.hpp>
#include <medusa/bits/approximations/Multiquadric.hpp>
#include <medusa/bits/approximations/InverseMultiquadric.hpp>

/**
 * @file
 * Instantiation of commonly used RBF bases.
 */

namespace mm {
template class RBFBasis<Gaussian<double>, Vec2d>;
template class RBFBasis<Multiquadric<double>, Vec2d>;
template class RBFBasis<InverseMultiquadric<double>, Vec2d>;
}  // namespace mm
