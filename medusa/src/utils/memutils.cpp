#include <medusa/bits/utils/memutils.hpp>

/**
 * @file
 * Implementation of memory related utilities.
 */

namespace mm {

std::string mem2str(std::size_t bytes) {
    double amount = bytes;
    std::vector<std::string> suffix = {"B", "kB", "MB", "GB"};
    for (int i = 0; i < 4; ++i) {
        if (amount < 100) {
            std::stringstream ss;
            ss << static_cast<int>(amount*10+0.5) / 10.0 << " " << suffix[i];
            return ss.str();
        }
        amount /= 1000;
    }
    return "More than your mem.";
}

}  // namespace mm
