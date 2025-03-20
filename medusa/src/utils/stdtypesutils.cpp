#include <medusa/bits/utils/stdtypesutils.hpp>

/**
 * @file
 * Implementation of utilities for std types.
 */

namespace mm {

std::vector<std::string> split(const std::string& str, const std::string& delim) {
    assert_msg(!delim.empty(), "Delimiter must not be empty.");
    std::vector<std::string> tokens;
    size_t prev = 0, pos = 0;
    do {
        pos = str.find(delim, prev);
        if (pos == std::string::npos) pos = str.length();
        tokens.emplace_back(str.substr(prev, pos-prev));
        prev = pos + delim.length();
    } while (pos < str.length() && prev < str.length());
    if (prev == str.length()) tokens.emplace_back("");
    return tokens;
}

std::vector<std::string> split(const std::string& str, char delim) {
    return split(str, std::string(1, delim));
}

/// @cond
std::string join(const std::vector<std::string>& parts, const std::string& joiner) {
    if (parts.empty()) return "";
    std::string result = parts[0];
    int n = parts.size();
    for (int i = 1; i < n; ++i) {
        result += joiner + parts[i];
    }
    return result;
}

std::string join(const std::vector<std::string>& parts, char joiner) {
    return join(parts, std::string(1, joiner));
}
/// @endcond

}  // namespace mm
