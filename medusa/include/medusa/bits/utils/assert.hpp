#ifndef MEDUSA_BITS_UTILS_ASSERT_HPP_
#define MEDUSA_BITS_UTILS_ASSERT_HPP_

/**
 * @file
 * Implementation of custom assert and debug utilities.
 */

#include <tinyformat/tinyformat.h>
#include "print.hpp"

namespace mm {

// print macro
/// @cond
#define VA_NUM_ARGS(...) VA_NUM_ARGS_IMPL(__VA_ARGS__, 5, 4, 3, 2, 1, 0)
#define VA_NUM_ARGS_IMPL(_1, _2, _3, _4, _5, N, ...) N
#define macro_dispatcher(func, ...)     macro_dispatcher_(func, VA_NUM_ARGS(__VA_ARGS__))
#define macro_dispatcher_(func, nargs)  macro_dispatcher__(func, nargs)
#define macro_dispatcher__(func, nargs) func ## nargs
#define addflag(a) {std::cerr << "flags=[flags, " << (a) << "];" << std::endl;}
#define prnv2(a, b) {std::cerr << a << " = " << (b) << ";" << std::endl;}
#define prnv1(a)   {std::cerr << #a << " = " << (a) << ";" << std::endl;}
/// @endcond
/**
 * Prints a variable name and value to standard output. Can take one or two parameters.
 * Example:
 * @code
 * int a = 6;
 * prn(a) // prints 'a = 6;'
 * prn("value", a) // prints 'value = 6;'
 * @endcode
 */
#define prn(...) macro_dispatcher(prnv, __VA_ARGS__)(__VA_ARGS__)

using tinyformat::printf;
using tinyformat::format;

/// Namespace holding custom assert implementation.
namespace assert_internal {
/**
 * Actual assert implementation.
 * @param condition Condition to test, e.g.\ `n > 0`.
 * @param file File where the assertion failed.
 * @param func_name Function name where the assertion failed.
 * @param line Line on which the assertion failed.
 * @param message Message as specified in the `assert_msg` macro.
 * @param format_list List of format field values to pass to `tinyformat::format` function.
 */
bool assert_handler_implementation(const char* condition, const char* file, const char* func_name,
                                   int line, const char* message, tfm::FormatListRef format_list);
/// Assert handler that unpacks varargs.
template<typename... Args>
bool assert_handler(const char* condition, const char* file, const char* func_name, int line,
                    const char* message, const Args&... args) {  // unpacks first argument
    tfm::FormatListRef arg_list = tfm::makeFormatList(args...);
    return assert_handler_implementation(condition, file, func_name, line, message, arg_list);
}
}  // namespace assert_internal

#ifdef NDEBUG
#define assert_msg(cond, ...) ((void)sizeof(cond))
#else
/**
 * @brief Assert with better error reporting.
 * @param cond Conditions to test.
 * @param ... The second parameter is also required and represents the message to print on failure.
 * For every %* field in message one additional parameter must be present.
 *
 * Example:
 * @code
 * assert_msg(n > 0, "n must be positive, got %d.", n);
 * @endcode
 */
#define assert_msg(cond, ...) ((void)(!(cond) && \
    mm::assert_internal::assert_handler( \
            #cond, __FILE__, __PRETTY_FUNCTION__, __LINE__, __VA_ARGS__) && (assert(0), 1)))
//            #cond, __FILE__, __PRETTY_FUNCTION__, __LINE__, __VA_ARGS__) && (exit(1), 1)))
#endif

/**
 * Prints given text in bold red.
 * @param s text to print.
 */
inline void print_red(const std::string& s) { std::cout << "\x1b[31;1m" << s << "\x1b[37;0m"; }
/**
 * Prints given text in bold white.
 * @param s text to print.
 */
inline void print_white(const std::string& s) { std::cout << "\x1b[37;1m" << s << "\x1b[37;0m"; }
/**
 * Prints given text in bold green.
 * @param s text to print.
 */
inline void print_green(const std::string& s) { std::cout << "\x1b[32;1m" << s << "\x1b[37;0m"; }

}  // namespace mm

#endif  // MEDUSA_BITS_UTILS_ASSERT_HPP_
