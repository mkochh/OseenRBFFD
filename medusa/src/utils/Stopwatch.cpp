#include <algorithm>
#include <medusa/bits/utils/assert.hpp>
#include <iomanip>
#include "medusa/bits/utils/Stopwatch.hpp"

/**
 * @file
 * Implementation of Stopwatch.
 */

namespace mm {

int Stopwatch::getID(const std::string& label) const {
    auto it = std::find(labels.begin(), labels.end(), label);
    assert_msg(it != labels.end(), "Label '%s' not found. Available labels: %s.",
               label.c_str(), labels);
    return it - labels.begin();
}
void Stopwatch::start(const std::string& label) {
    for (int i = 0; i < static_cast<int>(labels.size()); i++) {
        if (labels[i] == label) {
            assert_msg(!currently_running[i],
                       "start() for label '%s' was called more than once"
                       " before stop() call.", label.c_str());
            times[i] = std::chrono::steady_clock::now();
            currently_running[i] = true;
            return;
        }
    }
    // label not found
    labels.push_back(label);
    times.push_back(std::chrono::steady_clock::now());
    cumulative_time.push_back(0.0);
    counts.push_back(0);
    currently_running.push_back(true);
}
void Stopwatch::stop(const std::string& label) {
    int id = getID(label);
    assert_msg(currently_running[id],
               "stop() for label '%s' was called"
               " more than once before start() call.",
               label.c_str());
    ++counts[id];
    cumulative_time[id] += std::chrono::duration<double>(
            std::chrono::steady_clock::now()
            - times[getID(label)]).count();
    currently_running[id] = false;
}
double Stopwatch::cumulativeTime(const std::string& label) const {
    assert_msg(!isRunning(label),
               "Stopwatch with label '%s' is still running, `stop()` must be called"
               " before results can be displayed.", label.c_str());
    return cumulative_time[getID(label)];
}
int Stopwatch::numLaps(const std::string& label) const {
    assert_msg(!isRunning(label),
               "Stopwatch with label '%s' is still running, `stop()` must be called"
               " before results can be displayed.", label.c_str());
    return counts[getID(label)];
}
double Stopwatch::timePerLap(const std::string& label) const {
    assert_msg(!isRunning(label),
               "Stopwatch with label '%s' is still running, `stop()` must be called"
               " before results can be displayed.", label.c_str());
    return cumulativeTime(label)/static_cast<double>(numLaps(label));
}
void Stopwatch::clear() {
    times.clear();
    labels.clear();
    counts.clear();
    cumulative_time.clear();
    currently_running.clear();
}

/// Output average lap times for all labels.
std::ostream& operator<<(std::ostream& os, const Stopwatch& stopwatch) {
    size_t maxs = 0;
    for (const auto& label : stopwatch.labels) {
        if (label.size() > maxs) {
            maxs = label.size();
        }
    }
    for (const auto& label : stopwatch.labels) {
        os << std::setw(maxs) << label << ": " << stopwatch.timePerLap(label) << " [s]\n";
    }
    return os;
}

}  // namespace mm
