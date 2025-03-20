#ifndef MEDUSA_BITS_IO_CSV_FWD_HPP_
#define MEDUSA_BITS_IO_CSV_FWD_HPP_

/**
 * @file
 * Declarations for CSV I/O. To enable support for Eigen types, additional header CSV_Eigen.hpp
 * must be included.
 *
 * @example test/io/CSV_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/utils/assert.hpp>

/// @cond
namespace Eigen {
template <typename Derived>
class MatrixBase;

template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
class Matrix;
}
/// @endcond

namespace mm {

/**
 * Implements support for basic CSV I/O.
 *
 * @note To enable support for Eigen types, additional header CSV_Eigen.hpp must be included.
 *
 * @snippet io/CSV_test.cpp CSV usage example
 * @ingroup io
 */
class CSV {
  public:
    /**
     * Reads a CSV file to a vector of doubles. The file must have one number per line,
     * otherwise use CSV::read2d.
     * @param filename Path to the CSV file.
     * @throw Assertion fails if file cannot be read or if it is malformed.
     */
    static std::vector<double> read(const std::string& filename);

    /**
     * Reads a CSV file to a vector of vectors of doubles.
     * @param filename Path to the CSV file.
     * @param separator Character used as a separator in the CSV file.
     * @throw Assertion fails if file cannot be read or if it is malformed.
     */
    static std::vector<std::vector<double>> read2d(const std::string& filename,
                                                   char separator = ',');

    /**
     * Reads a CSV file to an Eigen matrix.
     * @param filename Path to the CSV file.
     * @param separator Character used as a separator in the CSV file.
     * @throw Assertion fails if file cannot be read or if it is malformed.
     * @note To use this function, additional header CSV_Eigen.hpp needs to be included.
     */
    static Eigen::Matrix<double, -1, -1, 0, -1, -1> readEigen(
            const std::string& filename, char separator = ',');

    /**
     * Writes given array to a CSV file (as a column).
     * @param filename Path to the CSV file.
     * @param array One dimensional array of data that supports range for loop iteration.
     * @throw Assertion fails if file cannot be written to.
     */
    template <typename arr_t>
    static void write(const std::string& filename, const arr_t& array);

    /**
     * Writes given 2d array to a CSV file.
     * @param filename Path to the CSV file.
     * @param array Two dimensional array of data that supports two nested range for loop
     * iterations.
     * @param separator Character to use as a separator in a CSV file.
     * @throw Assertion fails if file cannot be written to.
     */
    template <typename arr_t>
    static void write2d(const std::string& filename, const arr_t& array, char separator = ',');


    /**
     * Writes given Eigen matrix to a CSV file.
     * @param filename Path to the CSV file.
     * @param expr An Eigen expression that can be evaluated to a matrix.
     * @param separator Character to use as a separator in a CSV file.
     * @throw Assertion fails if file cannot be written to.
     */
    template <typename Derived>
    static void writeEigen(const std::string& filename, const Eigen::MatrixBase<Derived>& expr,
                           char separator = ',');

  private:
    /**
     * Splits given string into a vector of strings on the given separator.
     * Example:
     * @code
     * split("abc,def,efg", ',');  // returns {"abc", "def", "efg"}
     * @endcode
     */
    static std::vector<std::string> split(const std::string& str, char separator);
};

}  // namespace mm
#endif  // MEDUSA_BITS_IO_CSV_FWD_HPP_
