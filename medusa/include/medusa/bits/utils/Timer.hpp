#ifndef MEDUSA_BITS_UTILS_TIMER_HPP_
#define MEDUSA_BITS_UTILS_TIMER_HPP_

/**
 * @file
 * Timer class declarations.
 *
 * @example test/utils/Timer_test.cpp
 */

#include <medusa/Config.hpp>
#include <iostream>
#include <chrono>
#include <string>
#include <vector>

namespace mm {

/**
 * Simple timer class: add checkpoints throughout the code and measure
 * execution time between them.
 * When timing parts of code for example inside loops Stopwatch should be used instead.
 *
 * Usage example:
 * @snippet Timer_test.cpp Timer usage example
 * @sa Stopwatch
 * @ingroup utils
 */
class Timer {
  protected:
    typedef std::chrono::steady_clock::time_point time_type;  ///< Time type.
    std::vector<std::string> labels_;  ///< List of checkpoint labels.
    std::vector<time_type> times_;  ///< List of checkpoint times.

    /// Output timings between the checkpoints with given ids to `os`.
    void showTimings(int from, int to, std::ostream& os = std::cout) const;
    /// Returns the ID of checkpoints with a given label.
    int getID(const std::string& label) const;
    /// Return absolute time for a given `id`.
    time_type timeAt(int id) const;
  public:
    /// Adds a checkpoint with given label and remembers the time at which it was added.
    void addCheckPoint(const std::string& label);

    /**
     * Pretty print all durations between checkpoints. If there are less than 2 checkpoints,
     * only a warning is printed.
     */
    void showTimings(std::ostream& os = std::cout) const;
    /**
     * Output duration between the checkpoints with given labels to `os`.
     * Example: `showTimings("begin", "end");`
     */
    void showTimings(const std::string& from, const std::string& to,
                     std::ostream& os = std::cout) const;
    /// Return absolute time for a given label.
    time_type timeAt(const std::string& label) const;
    /// Return time difference in seconds between two checkpoints.
    double duration(const std::string& from, const std::string& to) const;
    /// Return time difference in seconds between `from` and now.
    double durationToNow(const std::string& from) const;
    /// Return the number of measurements taken.
    int size() const { return labels_.size(); }
    /// Return all labels.
    const std::vector<std::string>& labels() const { return labels_; }
    /// Return all times.
    const std::vector<time_type>& times() const { return times_; }
    /// Clear all of time points.
    void clear();

    /// Output all times between all checkpoints.
    friend std::ostream& operator<<(std::ostream& os, const Timer& timer) {
        timer.showTimings(os);
        return os;
    }
};

}  // namespace mm

#endif  // MEDUSA_BITS_UTILS_TIMER_HPP_
