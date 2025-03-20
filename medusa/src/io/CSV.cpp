#include "medusa/bits/io/CSV.hpp"

/**
 * @file
 * Implementation of CSV I/O utilities.
 */

namespace mm {

std::vector<double> CSV::read(const std::string& filename) {
    std::ifstream f(filename);
    assert_msg(f.good(), "Error opening CSV file '%s': %s.", filename, strerror(errno));

    std::vector<double> data;
    std::string s;
    while (getline(f, s)) {
        try {
            data.push_back(std::stod(s));
        } catch(const std::invalid_argument& e) {
            assert_msg(false, "Failed parsing '%s' to double.", s);
        }
    }
    return data;
}

std::vector<std::vector<double>> CSV::read2d(const std::string& filename, char separator) {
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

    std::vector<std::vector<double>> data(m, std::vector<double>(n));
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            try {
                data[i][j] = std::stod(lines[i][j]);
            } catch(const std::invalid_argument& e) {
                assert_msg(false, "Failed parsing '%s' to double.", lines[i][j]);
            }
        }
    }
    return data;
}

std::vector<std::string> CSV::split(const std::string& str, char separator) {
    std::vector<std::string> tokens;
    size_t prev = 0, pos = 0;
    do {
        pos = str.find(separator, prev);
        if (pos == std::string::npos) pos = str.length();
        tokens.emplace_back(str.substr(prev, pos-prev));
        prev = pos + 1;
    } while (pos < str.length() && prev < str.length());
    if (prev == str.length()) tokens.emplace_back("");
    return tokens;
}

}  // namespace mm
