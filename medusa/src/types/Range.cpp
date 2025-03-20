#include <medusa/bits/types/Range.hpp>
#include <medusa/bits/types/Vec.hpp>

/**
 * @file
 * Instantiation of common Range instances.
 */

template class mm::Range<double>;
template class mm::Range<int>;
/// @cond
template mm::Range<int> mm::Range<int>::seq(int);
template mm::Range<int> mm::Range<int>::seq(int, int);
template mm::Range<int> mm::Range<int>::seq(int, int, int);
/// @endcond
// template class mm::Range<bool>;
template class mm::Range<mm::Range<int>>;
template class mm::Range<mm::Vec1d>;
template class mm::Range<mm::Vec2d>;
template class mm::Range<mm::Vec3d>;
