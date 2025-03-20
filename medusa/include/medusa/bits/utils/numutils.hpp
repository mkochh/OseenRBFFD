#ifndef MEDUSA_BITS_UTILS_NUMUTILS_HPP_
#define MEDUSA_BITS_UTILS_NUMUTILS_HPP_

/**
 * @file
 * Declaration of numerical utilities.
 *
 * @example test/utils/numutils_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/types/Vec.hpp>
#include <medusa/bits/types/Range.hpp>
#include <cmath>

namespace mm {

/**
 * Ceils a floating point to an integer.
 * @warning This function might overflow if `x` is too large or return undefined results if
 * `x` is not a finite number.
 *
 * @ingroup utils
 */
template<typename T>
int iceil(T x) { return static_cast<int>(std::ceil(x)); }

/**
 * Floors a floating point to an integer.
 * @warning This function might overflow if `x` is too large or return undefined results if
 * `x` is not a finite number.
 *  @ingroup utils
 */
template<typename T>
int ifloor(T x) { return static_cast<int>(std::floor(x)); }

/** Compile time integer power, returns `base` raised to power `exponent`.  @ingroup utils */
template <unsigned int exponent>
inline double ipow(double base) {
    return ipow<exponent - 1>(base) * base;
}

/** Compile time integer power (base case 0)  @ingroup utils */
template <>
inline double ipow<0>(double) { return 1; }

/**
 * Compute non-negative integer power `b^e`. This function is usually faster
 * than `std::pow` for `e < 50` and matches the `std::pow` performance at
 * approximately `e = 100`. For compile time constant `e` this function
 * use @ref ipow.
 * @ingroup utils
 */
template <typename T>
inline T ipow(T b, int e) {
    T r = 1.0;
    while (e > 0) {
        r *= b;
        --e;
    }
    return r;
}

/**
 * Compute possibly negative integer power `b^e`. This function is usually faster
 * than `std::pow` for `e < 50` and matches the `std::pow` performance at
 * approximately `e = 100`.
 * @ingroup utils
 */
template <typename T>
inline T ipowneg(T b, int e) {
    if (e < 0) return 1.0 / ipow(b, -e);
    else return ipow(b, e);
}

/// Signum overload for unsigned types
template <typename T>
inline constexpr int signum(T x, std::false_type) {
    return T(0) < x;
}
/// Signum overload for unsigned types
template <typename T>
inline constexpr int signum(T x, std::true_type) {
    return (T(0) < x) - (x < T(0));
}
/**
 * Signum function: determines a sign of a number `x`.
 * @param x A number under inspection.
 * @return `0` if `x == 0`, `-1` if `x` is negative and `+1` if `x` is positive
 * @ingroup utils
 */
template <typename T>
inline constexpr int signum(T x) {
    return signum(x, std::is_signed<T>());
}

/**
 * Increments a multi-dimensional  counter with given limits.
 * @param counter A valid counter state.
 * @param limit Limit for each dimension of a counter.
 * @return `true` if an increment was performed and `false` otherwise.
 * Repeated application of this function to counter with initial state
 * `0 0 0` and limits `1 2 3` yields a sequence:
 * @code
 * 0 0 0 -> true
 * 0 0 1 -> true
 * 0 0 2 -> true
 * 0 1 0 -> true
 * 0 1 1 -> true
 * 0 1 2 -> true
 * 0 1 2 -> false
 * @endcode
 * This is made so that it can be used nicely with a `do {} while ()` loop, similarly to
 * `std::next_permutation`.
 * @snippet test/utils/numutils_test.cpp Increment counter usage
 * @ingroup utils
 */
template <int dim>
bool incrementCounter(Vec<int, dim>& counter, const Vec<int, dim>& limit) {
    for (int i = dim - 1; i >= 0; --i) {
        if (counter[i] >= limit[i] - 1) {
            counter[i] = 0;
        } else {
            counter[i]++;
            return true;
        }
    }
    return false;
}

/**
 * Increments a multi-dimensional counter with given upper and lower limits.
 * @param counter A valid counter state, which is modified.
 * @param low Lower limit for each dimension of a counter.
 * @param high Upper limit for each dimension of a counter. It must hold that
 * `low <= high` in each component.
 * @return `true` if an increment was performed and `false` otherwise.
 * @sa incrementCounter
 * @ingroup utils
 */
template <int dim>
bool incrementCounter(Vec<int, dim>& counter, const Vec<int, dim>& low, const Vec<int, dim>& high) {
    for (int i = dim - 1; i >= 0; --i) {
        if (counter[i] >= high[i] - 1) {
            counter[i] = low[i];
        } else {
            counter[i]++;
            return true;
        }
    }
    return false;
}

/**
 * Increments a multi-dimensional counter with given upper and lower limits and global upper size,
 * looping around to 0 if needed. It need <b>not</b> hold that `low <= high` in each component.
 * @param counter A valid counter state, which is modified.
 * @param low Lower limit for each dimension of a counter.
 * @param high Upper limit for each dimension of a counter.
 * @param size Global upper limit.
 * @return `true` if an increment was performed and `false` otherwise.
 * Repeated application of this function to counter with limits
 * `1 0 2` and `0 2 1` with size limit `3 3 3`, starting from lower limit yields a sequence:
 * @code
 * 1 0 2 -> true
 * 1 0 0 -> true
 * 1 1 2 -> true
 * 1 1 0 -> true
 * 2 0 2 -> true
 * 2 0 0 -> true
 * 2 1 2 -> true
 * 2 1 0 -> true
 * 2 1 0 -> false
 * @endcode
 * This is useful for cyclic iteration over multidimensional arrays.
 * @snippet test/utils/numutils_test.cpp Cyclic counter usage
 * @sa incrementCounter
 * @ingroup utils
 */
template <int dim>
bool incrementCyclicCounter(Vec<int, dim>& counter, const Vec<int, dim>& low,
                            const Vec<int, dim>& high, const Vec<int, dim>& size) {
    for (int i = dim - 1; i >= 0; --i) {
        if (counter[i] == high[i] - 1 || (high[i] == 0 && counter[i] == size[i] - 1)) {
            counter[i] = low[i];
        } else {
            counter[i]++;
            if (counter[i] >= size[i]) counter[i] = 0;
            return true;
        }
    }
    return false;
}

/**
 * Multidimensional clone of Matlab's linspace function.
 * Uniformly discretizes cuboid given with `beg` and `end` points.
 * Similar to numpy's linspace and Matlab's meshgrid.
 * @param beg Beginning of a cuboid.
 * @param end Ending of a cuboid.
 * @param counts How many discretization points to use in each dimension.
 * @param include_boundary Flag whether to include boundary of a cuboid in a given
 * dimension.
 * @return Uniform discretization points of a domain with counts[i] points in dimension i.
 * Instead of returning a n-dim matrix it returns it as a Range instead.
 * @ingroup utils
 */
template <class scalar_t, int dim>
Range<Vec<scalar_t, dim>> linspace(const Vec<scalar_t, dim>& beg, const Vec<scalar_t, dim>& end,
                                   const Vec<int, dim>& counts,
                                   const Vec<bool, dim> include_boundary = true) {
    Range<Vec<scalar_t, dim>> ret;
    for (int i = 0; i < dim; ++i) {
        if (include_boundary[i]) assert(counts[i] >= 2);
        else assert(counts[i] >= 0);
        if (counts[i] == 0) return ret;
    }
    Vec<scalar_t, dim> step;
    for (int i = 0; i < dim; ++i) {
        if (include_boundary[i]) step[i] = (end[i] - beg[i]) / (counts[i] - 1);
        else step[i] = (end[i] - beg[i]) / (counts[i] + 1);
    }
    Vec<int, dim> counter = 0;
    do {
        Vec<scalar_t, dim> tmp;
        for (int i = 0; i < dim; ++i) {
            if (include_boundary[i]) tmp[i] = beg[i] + counter[i] * step[i];
            else tmp[i] = beg[i] + (counter[i] + 1) * step[i];
        }
        ret.push_back(tmp);
    } while (incrementCounter(counter, counts));
    return ret;
}

/// Overload for bool argument of `include_boundary`. @ingroup utils
template <class scalar_t, int dim>
Range<Vec<scalar_t, dim>> linspace(const Vec<scalar_t, dim>& beg, const Vec<scalar_t, dim>& end,
                                   const Vec<int, dim>& counts, bool include_boundary) {
    return linspace(beg, end, counts, Vec<bool, dim>(include_boundary));
}

/// Overload for 1 dimension. @ingroup utils
template <typename scalar_t>
Range<scalar_t> linspace(scalar_t beg, scalar_t end, int count, bool include_boundary = true) {
    Range<scalar_t> ret;
    if (include_boundary) {
        assert_msg(count >= 2, "Count must be >= 2, got %d.", count);
        scalar_t step = (end - beg) / (count - 1);
        for (int i = 0; i < count; ++i) {
            ret.push_back(beg + step * i);
        }
    } else {
        assert_msg(count >= 0, "Count must be >= 0, got %d.", count);
        scalar_t step = (end - beg) / (count + 1);
        for (int i = 1; i <= count; ++i) {
            ret.push_back(beg + step * i);
        }
    }
    return ret;
}

/**
 * Solves `f(x) = target` using bisection.
 * @tparam function_t Function type, like std::function or lambda, that returns output_type.
 * @tparam input_t Floating point data type that support arithmetic operations.
 * @tparam output_t Floating point data type, such that `mm::signum(output_t)` can be called.
 * @tparam verbose Reports approximation every iteration.
 * @param f Function mapping `input_t -> output_t`, for which to solve `f(x) = target`.
 * @param lo Lower bound of the interval, containing the solution.
 * @param hi Upper bound of the interval, containing the solution.
 * @param target Target value, default 0.
 * @param tolerance Desired accuracy of the solution, default `1e-4`.
 * @param max_iter Maximal number of iterations.
 * @return The solution `x`, such that `f(x)` is approximately equal to `target`.
 * @ingroup utils
 */
template<typename function_t, typename input_t, typename output_t, bool verbose = false>
input_t bisection(const function_t& f, input_t lo, input_t hi, output_t target = 0.0,
                  input_t tolerance = 1e-4, int max_iter = 40) {
    input_t a = lo, d = hi - lo;
    output_t fa = f(lo) - target;
    output_t fb = f(hi) - target;

    assert_msg(signum(fa) != signum(fb),
               "Function values have the same sign, f(a) = %s, f(b) = %s.", fa, fb);

    while (d > tolerance) {
        d /= static_cast<input_t>(2.);
        output_t fc = f(a+d) - target;

        if (verbose) std::cout << "f(" << a+d << ")=" << fc + target << "\n";
        if (signum(fa) == signum(fc)) {
            a += d;
            fa = fc;
        }
        --max_iter;
        if (max_iter < 0) {
            std::cerr << "\nBisection ended without finding the solution after max_num"
                    " iterations." << std::endl;
            break;
        }
    }
    return a+d;
}

}  // namespace mm

#endif  // MEDUSA_BITS_UTILS_NUMUTILS_HPP_
