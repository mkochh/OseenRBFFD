#include <medusa/bits/spatial_search/Grid.hpp>

/**
 * @file
 * Instantiations of Grid.
 */

/// @cond
template class mm::Grid<int, 1>;
template std::ostream& mm::operator<<(std::ostream&, const mm::Grid<int, 1>&);
template class mm::Grid<int, 2>;
template std::ostream& mm::operator<<(std::ostream&, const mm::Grid<int, 2>&);
template class mm::Grid<int, 3>;
template std::ostream& mm::operator<<(std::ostream&, const mm::Grid<int, 3>&);
/// @endcond
