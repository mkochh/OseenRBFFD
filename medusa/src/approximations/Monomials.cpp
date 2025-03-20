#include <medusa/bits/approximations/Monomials.hpp>
#include <medusa/bits/types/Vec.hpp>

/**
 * @file
 * Instantiations of most commonly used Monomial bases.
 */

/// @cond
template class mm::Monomials<mm::Vec1d>;
template class mm::Monomials<mm::Vec2d>;
template class mm::Monomials<mm::Vec3d>;
template std::ostream& mm::operator<<(std::ostream& os, const mm::Monomials<mm::Vec1d>&);
template std::ostream& mm::operator<<(std::ostream& os, const mm::Monomials<mm::Vec2d>&);
template std::ostream& mm::operator<<(std::ostream& os, const mm::Monomials<mm::Vec3d>&);
/// @endcond
