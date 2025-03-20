#include "medusa/bits/approximations/Polyharmonic.hpp"

/**
 * @file
 * Explicit instantiations for common parameters.
*/

/// @cond
template class mm::Polyharmonic<double>;
template double mm::Polyharmonic<double>::operator()(double r2, mm::Lap<1>) const;
template double mm::Polyharmonic<double>::operator()(double r2, mm::Lap<2>) const;
template double mm::Polyharmonic<double>::operator()(double r2, mm::Lap<3>) const;
template std::ostream& mm::operator<<(std::ostream& os, const mm::Polyharmonic<double>&);

template class mm::Polyharmonic<double, 3>;
template double mm::Polyharmonic<double, 3>::operator()(double r2, mm::Lap<1>) const;
template double mm::Polyharmonic<double, 3>::operator()(double r2, mm::Lap<2>) const;
template double mm::Polyharmonic<double, 3>::operator()(double r2, mm::Lap<3>) const;
template std::ostream& mm::operator<<(std::ostream& os, const mm::Polyharmonic<double, 3>&);

template class mm::Polyharmonic<double, 5>;
template double mm::Polyharmonic<double, 5>::operator()(double r2, mm::Lap<1>) const;
template double mm::Polyharmonic<double, 5>::operator()(double r2, mm::Lap<2>) const;
template double mm::Polyharmonic<double, 5>::operator()(double r2, mm::Lap<3>) const;
template std::ostream& mm::operator<<(std::ostream& os, const mm::Polyharmonic<double, 5>&);

template class mm::Polyharmonic<double, 7>;
template double mm::Polyharmonic<double, 7>::operator()(double r2, mm::Lap<1>) const;
template double mm::Polyharmonic<double, 7>::operator()(double r2, mm::Lap<2>) const;
template double mm::Polyharmonic<double, 7>::operator()(double r2, mm::Lap<3>) const;
template std::ostream& mm::operator<<(std::ostream& os, const mm::Polyharmonic<double, 7>&);
/// @endcond
