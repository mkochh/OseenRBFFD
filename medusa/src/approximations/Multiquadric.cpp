#include "medusa/bits/approximations/Multiquadric.hpp"

/**
 * @file
 * Explicit instantiations for common parameters.
 */

/// @cond
template class mm::Multiquadric<double>;
template double mm::Multiquadric<double>::operator()(double r2, mm::Lap<1>) const;
template double mm::Multiquadric<double>::operator()(double r2, mm::Lap<2>) const;
template double mm::Multiquadric<double>::operator()(double r2, mm::Lap<3>) const;
template std::ostream& mm::operator<<(std::ostream& os, const mm::Multiquadric<double>&);
/// @endcond
