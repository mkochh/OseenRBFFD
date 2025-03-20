#include "medusa/bits/utils/Stopwatch.hpp"
#include <thread>

#include "gtest/gtest.h"

namespace mm {

TEST(Utils, Stopwatch) {
    double max_timer_error = 1e-4;
    Stopwatch t;
    // Same as Timer.duration test.
    std::chrono::milliseconds sleep_time(1);
    std::chrono::steady_clock::time_point t1, t2;
    t1 = std::chrono::steady_clock::now();
    t.start("a");
    t.stop("a");
    t2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> time_span =
            std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    double diff = t.cumulativeTime("a");
    EXPECT_GT(diff, 0);
    EXPECT_LT(diff, time_span.count());
    EXPECT_EQ(t.numLaps("a"), 1);

    // Testing loop.
    int n_loops = 10;
    t1 = std::chrono::steady_clock::now();
    for (int i = 0; i < n_loops; i++) {
        t.start("b");
        std::this_thread::sleep_for(sleep_time);
        t.stop("b");
    }
    t2 = std::chrono::steady_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    double per_lap_ref = time_span.count() / n_loops;
    EXPECT_EQ(n_loops, t.numLaps("b"));
    double per_lap = t.timePerLap("b");
    EXPECT_LT(per_lap, per_lap_ref);
    EXPECT_LT(per_lap_ref - per_lap, max_timer_error);
    // Testing loop - two stopwatches. The enclosing one is assumed
    // to be correct since the previous test did pass.
    for (int i = 0; i < n_loops; i++) {
        t.start("d");
        t.start("c");
        std::this_thread::sleep_for(sleep_time);
        t.stop("c");
        t.start("c");
        std::this_thread::sleep_for(sleep_time);
        t.stop("c");
        t.stop("d");
    }
    EXPECT_EQ(2 * n_loops, t.numLaps("c"));
    double per_lap_ref_2 = t.timePerLap("d") / 2;
    double per_lap_2 = t.timePerLap("c");
    EXPECT_LT(per_lap_2, per_lap_ref_2);
    EXPECT_LT(per_lap_ref_2 - per_lap_2, 2 * max_timer_error);
    t.clear();
}

TEST(Utils, DISABLED_StopwatchUsageExample) {
    /// [Stopwatch usage example]
    Stopwatch s;
    for (int i = 0; i < 10; ++i) {
        // ... code: NOT timed  ...
        s.start("section1");
        // ... code: timed  ...
        s.stop("section1");
        // ... code: NOT timed  ...
        s.start("section2");
        // ... code: timed  ...
        s.stop("section2");
    }
    s.cumulativeTime("section1");  // Returns cumulative time.
    s.numLaps("section1");  // Returns number of laps counted.
    s.timePerLap("section1");  // Returns average time per lap.

    std::cout << s << std::endl;
    /// [Stopwatch usage example]
}

}  // namespace mm
