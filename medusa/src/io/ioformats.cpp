#include <medusa/bits/io/ioformat.hpp>

/**
 * @file
 * Declarations of most common formats.
 */

namespace mm {
/// @cond
Eigen::IOFormat CleanFmt(4, 0, ", ", ",\n", "[", "]");
Eigen::IOFormat InlineFmt(Eigen::StreamPrecision, Eigen::DontAlignCols,
                          ", ", "; ", "[", "]", "", "");
Eigen::IOFormat MatlabFmt(Eigen::FullPrecision, Eigen::DontAlignCols,
                          ", ", "; ", "", "", "[", "];");
Eigen::IOFormat MathematicaFmt(Eigen::FullPrecision, Eigen::DontAlignCols,
                               ", ", ", ", "{", "}", "{", "};");
Eigen::IOFormat CSVFmt(Eigen::FullPrecision, 0, ", ");
/// @endcond
}  // namespace mm
