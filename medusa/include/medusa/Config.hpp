#ifndef MEDUSA_CONFIG_HPP_
#define MEDUSA_CONFIG_HPP_

/**
 * @file
 * Contains compile time settings for the project, such as Eigen plugins.
 */

/// Set Eigen MatrixBase plugin
#define EIGEN_MATRIXBASE_PLUGIN "medusa/bits/types/MatrixBaseAddons.hpp"
/// Set Eigen Matrix plugin
#define EIGEN_MATRIX_PLUGIN "medusa/bits/types/MatrixAddons.hpp"
/// Set Eigen default IO format

#ifdef EIGEN_DEFAULT_IO_FORMAT
#error "Medusa's Eigen configuration was included after Eigen. Make sure you include this header "
       "before including any of the Eigen headers."
# include <stophereandnow>  // Make sure to not compile further and emit a bunch of Eigen errors.
#else
#define EIGEN_DEFAULT_IO_FORMAT Eigen::IOFormat(Eigen::StreamPrecision, \
                                                Eigen::DontAlignCols, ", ", "; ", "", "", "[", "]")
#endif
/// Strong inline macro
#define SINL EIGEN_STRONG_INLINE

#ifndef NDEBUG
/// Initialize matrices with NaN in debug mode for easier error detection.
#define EIGEN_INITIALIZE_MATRICES_BY_NAN
#endif

#include <vector>

/// Root namespace for the whole library.
namespace mm {

typedef std::vector<int> indexes_t;  ///< Class representing a collection of indices.

/// Value of Pi in type `T`. Usage: @code double tau = 2*Pi<double>::value;@endcode @ingroup utils
template <typename T>
struct Pi {
    static constexpr T value = T(3.14159265358979323846264338327950L);  ///< Value of @f$\pi@f$.
};

static const double PI = Pi<double>::value;  ///< Mathematical constant pi in double precision.
static const double INF = 1.0 / 0.0;  ///< Infinite floating point value.
static const double NaN = 0.0 / 0.0;  ///< Not-a-number floating point value.

}  // namespace mm

#endif  // MEDUSA_CONFIG_HPP_
