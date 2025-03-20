#ifndef MEDUSA_BITS_IO_IOFORMAT_HPP_
#define MEDUSA_BITS_IO_IOFORMAT_HPP_

/**
 * @file
 * File defining output formats for Eigen matrices.
 */

#include <Eigen/Core>

namespace mm {

/// Clean readable multiline aligned format.
extern Eigen::IOFormat CleanFmt;
/// Readable inline format.
extern Eigen::IOFormat InlineFmt;
/// Full precision format understood by Matlab/Octave.
extern Eigen::IOFormat MatlabFmt;
/// Full precision output understood by Wolfram Mathematica.
extern Eigen::IOFormat MathematicaFmt;
/// Valid CSV format.
extern Eigen::IOFormat CSVFmt;

}  // namespace mm

#endif  // MEDUSA_BITS_IO_IOFORMAT_HPP_
