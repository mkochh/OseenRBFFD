#ifndef MEDUSA_BITS_UTILS_RANDUTILS_HPP_
#define MEDUSA_BITS_UTILS_RANDUTILS_HPP_

/**
 * @file
 * Utilities for randomization.
 *
 * @example test/utils/randutils_test.cpp
 */

#include <medusa/Config.hpp>
#include <medusa/bits/utils/assert.hpp>
#include <medusa/bits/types/Range.hpp>
#include <random>
#include <chrono>

namespace mm {

/**
 * Return a random seed. The seed is truly random if available, otherwise it is
 * the current system time.
 *
 * @ingroup utils
 */
inline unsigned int get_seed() {
    try {
        std::random_device rd;
        return rd();
    } catch (std::exception& e) {
        return std::chrono::system_clock::now().time_since_epoch().count();
    }
}

template <typename T, typename URNG>
T random_choice(const Range<T>& elements, const Range<double>& weights, bool normed,
                URNG& generator);

/**
 * Randomly returns one of the specified elements with distribution according to given
 * weights. A `std::mt19937` generator is created and seeded randomly with call to get_seed().
 * @param elements A pool of elements to choose from.
 * @param weights Weights of the elements. If this argument is omitted all elements are
 * assigned the same weights.
 * @param normed Boolean indicating that the weights are already normed and that
 * additional computation is not necessary.
 * @throws Assertion might fail if thw weights are claimed to be normalized but are not.
 * @sa get_seed
 * @ingroup utils
 */
template <typename T>
T random_choice(const Range<T>& elements, const Range<double>& weights = {},
                bool normed = false) {
    std::mt19937 generator(get_seed());
    return random_choice(elements, weights, normed, generator);
}

/// Overload for custom generator. @sa random_choice @ingroup utils
template <typename T, typename URNG>
T random_choice(const Range<T>& elements, const Range<double>& weights, bool normed,
                URNG& generator) {
    Range<double> actual_weights = weights;
    if (actual_weights.empty()) {
        normed = true;
        actual_weights.assign(elements.size(), 1.0 / elements.size());
    }
    assert_msg(actual_weights.size() == elements.size(),
               "Weights not specified for all elements. Got %d weights but %d elements.",
               weights.size(), elements.size());
    if (!normed) {  // must be 0.0 so that the inferred type is double
        double sum = std::accumulate(actual_weights.begin(), actual_weights.end(), 0.0);
        for (auto& w : actual_weights) w /= sum;
    }
    std::uniform_real_distribution<double> uniform(0, 1);
    double x = uniform(generator);
    for (int i = 0; i < elements.size(); ++i) {
        if (x < actual_weights[i]) return elements[i];
        x -= actual_weights[i];
    }
    assert_msg(false, "Should never get here, you might have claimed the weights "
                      "are normalized, but weren't.");
    return elements[0];
}

}  // namespace mm

#endif  // MEDUSA_BITS_UTILS_RANDUTILS_HPP_
