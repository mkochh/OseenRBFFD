/**
 * @file
 * Timer class implementations.
 */

#include <medusa/bits/utils/Timer.hpp>
#include <medusa/bits/utils/assert.hpp>
#include <iomanip>

namespace mm {

void Timer::addCheckPoint(const std::string& label) {
    assert_msg(std::find(labels_.begin(), labels_.end(), label) == labels_.end(),
               "Label '%s' already exists. Use stopwatch to time repeatedly.", label);
    labels_.push_back(label);
    times_.push_back(std::chrono::steady_clock::now());
}
void Timer::showTimings(std::ostream& os) const {
    if (labels_.size() < 2) {
        os << "Not enough checkpoints.";
        return;
    }
    size_t M = 0;
    for (const auto& label : labels_) M = std::max(M, label.size());
    for (size_t c = 1; c < labels_.size(); ++c) {
        os << std::left
           << std::setw(M) << labels_[c - 1] << " -- " << std::setw(M) << labels_[c]
           << ' ' << std::setw(10) << std::scientific
           << std::chrono::duration<double>(times_[c] - times_[c-1]).count() << " [s]"
           << std::endl;
    }
    os << std::left << std::setw(2*M+5)
       << "total time " << std::setw(10)  // << std::scientific
       << std::chrono::duration<double>(times_.back() - times_[0]).count()
       << " [s]" << std::endl;
}

void Timer::showTimings(const std::string& from, const std::string& to,
                        std::ostream& os) const {
    showTimings(getID(from), getID(to), os);
}

void Timer::showTimings(int from, int to, std::ostream& os) const {
    assert_msg(0 <= from && from < size(), "From ID %d out of range [0, %d).", from, size());
    assert_msg(0 <= to && to < size(), "To ID %d out of range [0, %d).", to, size());
    os << std::left << std::setw(20) << labels_[from] << " -- "
       << std::setw(20) << labels_[to] << std::setw(18)
       << std::chrono::duration<double>(times_[to] - times_[from]).count()
       << "[s]" << std::endl;
}

Timer::time_type Timer::timeAt(int id) const {
    assert_msg(0 <= id && id < size(), "ID %d out of range [0, %d).", id, size());
    return times_[id];
}

Timer::time_type Timer::timeAt(const std::string& label) const {
    return times_[getID(label)];
}
double Timer::duration(const std::string& from, const std::string& to) const {
    return std::chrono::duration<double>(times_[getID(to)] - times_[getID(from)]).count();
}
double Timer::durationToNow(const std::string& from) const {
    return std::chrono::duration<double>(
            std::chrono::steady_clock::now() - times_[getID(from)]).count();
}
int Timer::getID(const std::string& label) const {
    auto it = std::find(labels_.begin(), labels_.end(), label);
    assert_msg(it != labels_.end(), "Label '%s' not found. Available labels: %s.",
               label.c_str(), labels_);
    return it - labels_.begin();
}
void Timer::clear() {
    times_.clear();
    labels_.clear();
}

}  // namespace mm
