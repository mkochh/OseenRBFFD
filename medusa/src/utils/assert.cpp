#include <medusa/bits/utils/assert.hpp>

/**
 * @file
 * Implementation of custom assert utilities.
 */

namespace mm {

namespace assert_internal {

bool assert_handler_implementation(const char* condition, const char* file, const char* func_name,
                                   int line, const char* message, tfm::FormatListRef format_list) {
    std::cerr << "\x1b[37;1m";  // white bold
    tfm::format(std::cerr, "%s:%d: %s: Assertion `%s' failed with message:\n",
                file, line, func_name, condition);
    std::cerr << "\x1b[31;1m";  // red bold
    tfm::vformat(std::cerr, message, format_list);
    std::cerr << "\x1b[37;0m\n";  // no color
    std::cerr.flush();
    return true;
}

}  // namespace assert_internal

}  // namespace mm
