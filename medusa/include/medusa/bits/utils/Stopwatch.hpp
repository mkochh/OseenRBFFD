#ifndef MEDUSA_BITS_UTILS_STOPWATCH_HPP_
#define MEDUSA_BITS_UTILS_STOPWATCH_HPP_

/**
 * @file
 * Declaration of Stopwatch.
 *
 * @example test/utils/Stopwatch_test.cpp
 */

#include <medusa/Config.hpp>
#include <vector>
#include <chrono>
#include <string>
#include <ostream>

namespace mm {

/**
 * A simple stopwatch class: time sections of code that execute repeatedly and get average
 * execution time. Stopwatch keeps track of cumulative time spent between all calls to
 * `start(label)`, `stop(label)` pairs for the same `label`, number of laps
 * (calls of start/stop pairs), and cumulative time per lap.
 *
 * For timing execution time between various checkpoints, use the Timer class.
 *
 * @snippet Stopwatch_test.cpp Stopwatch usage example
 * @sa Timer
 * @ingroup utils
 */
class Stopwatch {
    typedef std::chrono::steady_clock::time_point time_type;  ///< Time type.
    std::vector<std::string> labels;  ///< List of stopwatch labels.
    std::vector<double> cumulative_time;  ///< List of cumulative times for each stopwatch.
    std::vector<int> counts;  ///< List of lap counts for each stopwatch.
    std::vector<time_type> times;  ///< List of times for each stopwatch.
    std::vector<bool> currently_running;  ///< For tracking when stopwatch is running.

    /// Returns the ID of stopwatch with a given label.
    int getID(const std::string& label) const;

  public:
    /// Starts the stopwatch with a given label.
    void start(const std::string& label);
    /**
     * Stops the stopwatch with a given label.
     * Updates cumulative time and number of laps information for the given label.
     * Individual lap durations are not saved.
     */
    void stop(const std::string& label);
    /// Returns total time of all laps for a given label.
    double cumulativeTime(const std::string& label) const;
    /// Returns number of laps for a given label.
    int numLaps(const std::string& label) const;
    /// Returns average time spent per lap.
    double timePerLap(const std::string& label) const;
    /// Returns if stopwatch with a given label is currently running.
    bool isRunning(const std::string& label) const { return currently_running[getID(label)]; }
    /// Clear all stopwatch related data.
    void clear();
    /// Output average lap times for all labels.
    friend std::ostream& operator<<(std::ostream& os, const Stopwatch& stopwatch);
};

}  // namespace mm

#endif  // MEDUSA_BITS_UTILS_STOPWATCH_HPP_
