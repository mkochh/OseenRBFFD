#include "medusa/bits/approximations/InverseMultiquadric.hpp"

/**
 * @file
 * Explicit instantiations for common parameters.
 */

/// @cond
template class mm::InverseMultiquadric<double>;
template double mm::InverseMultiquadric<double>::operator()(double r2, mm::Lap<1>) const;
template double mm::InverseMultiquadric<double>::operator()(double r2, mm::Lap<2>) const;
template double mm::InverseMultiquadric<double>::operator()(double r2, mm::Lap<3>) const;
template std::ostream& mm::operator<<(std::ostream& os, const mm::InverseMultiquadric<double>&);
/// @endcond
