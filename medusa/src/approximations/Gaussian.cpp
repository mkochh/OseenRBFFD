#include <medusa/bits/approximations/Gaussian.hpp>

/**
 * @file
 * Instantiation of most commonly used Gaussian RBF.
 */

/// @cond
template class mm::Gaussian<double>;
template double mm::Gaussian<double>::operator()(double r2, mm::Lap<1>) const;
template double mm::Gaussian<double>::operator()(double r2, mm::Lap<2>) const;
template double mm::Gaussian<double>::operator()(double r2, mm::Lap<3>) const;
template std::ostream& mm::operator<<(std::ostream& os, const mm::Gaussian<double>&);
/// @endcond
