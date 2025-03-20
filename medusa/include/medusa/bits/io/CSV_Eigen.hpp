#ifndef MEDUSA_BITS_IO_CSV_EIGEN_HPP_
#define MEDUSA_BITS_IO_CSV_EIGEN_HPP_

/**
 * @file
 * Implementation of CSV I/O utilities for Eigen types.
 */

#include "CSV_fwd.hpp"
#include <Eigen/Core>

namespace mm {

inline Eigen::MatrixXd CSV::readEigen(const std::string& filename, char separator) {
    std::ifstream f(filename);
    assert_msg(f.good(), "Error opening CSV file '%s': %s.", filename, strerror(errno));

    std::vector<std::vector<std::string>> lines;
    std::string s;
    while (getline(f, s)) lines.push_back(split(s, separator));

    size_t m = lines.size();
    size_t n = lines[0].size();
    for (size_t i = 1; i < m; ++i) {
        assert_msg(lines[i].size() == n, "Not all rows in CSV file have the same number of "
                                         "elements, row 1 has %d and row %d has %d.",
                   n, i, lines[i].size());
    }

    Eigen::MatrixXd M(m, n);
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            try {
                M(i, j) = std::stod(lines[i][j]);
            } catch(const std::invalid_argument& e) {
                assert_msg(false, "Failed parsing '%s' to double.", lines[i][j]);
            }
        }
    }
    return M;
}

template <typename Derived>
void CSV::writeEigen(const std::string& filename, const Eigen::MatrixBase <Derived>& expr,
                     char separator) {
    std::ofstream f(filename);
    assert_msg(f.good(), "Error opening CSV file '%s': %s.", filename, strerror(errno));
    f << std::setprecision(16);
    int r = expr.rows(), c = expr.cols();
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            if (j > 0) f << separator;
            f << expr(i, j);
        }
        f << '\n';
    }
}

}  // namespace mm

#endif  // MEDUSA_BITS_IO_CSV_EIGEN_HPP_
