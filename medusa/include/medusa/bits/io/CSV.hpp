#ifndef MEDUSA_BITS_IO_CSV_HPP_
#define MEDUSA_BITS_IO_CSV_HPP_

/**
 * @file
 * Implementation of CSV I/O utilities.
 */

#include "CSV_fwd.hpp"
#include <fstream>
#include <cstring>
#include <exception>
#include <string>
#include <iomanip>

namespace mm {

template <typename arr_t>
void CSV::write(const std::string& filename, const arr_t& array) {
    std::ofstream f(filename);
    assert_msg(f.good(), "Error opening CSV file '%s': %s.", filename, strerror(errno));
    f << std::setprecision(16);
    for (const auto& x : array) {
        f << x << '\n';
    }
}

template <typename arr_t>
void CSV::write2d(const std::string& filename, const arr_t& array, char separator) {
    std::ofstream f(filename);
    assert_msg(f.good(), "Error opening CSV file '%s': %s.", filename, strerror(errno));
    f << std::setprecision(16);
    for (const auto& x : array) {
        bool first = true;
        for (const auto& y : x) {
            if (!first) f << separator;
            else first = false;
            f << y;
        }
        f << '\n';
    }
}

}  // namespace mm

#endif  // MEDUSA_BITS_IO_CSV_HPP_
