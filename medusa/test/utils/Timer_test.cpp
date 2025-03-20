#include <medusa/bits/utils/Timer.hpp>

#include "gtest/gtest.h"

namespace mm {

TEST(Utils, TimerDuration) {
    Timer t;
    t.addCheckPoint("a");
    t.addCheckPoint("b");
    double diff = t.duration("a", "b");
    EXPECT_GT(diff, 0);
    double diff2 = t.durationToNow("a");
    EXPECT_LT(diff, diff2);
    EXPECT_GT(diff2, 0);
}

TEST(Utils, DISABLED_TimerUsageExample) {
    /// [Timer usage example]
    Timer timer;
    timer.addCheckPoint("beg");
    //... code ...
    timer.addCheckPoint("mid");
    // ...more code ...
    timer.addCheckPoint("end");

    timer.showTimings("beg", "mid");  // shows time between two checkpoints
    std::cout << timer << std::endl;

    double t = timer.duration("beg", "end");  // time difference in seconds
    /// [Timer usage example]
    std::cout << t << std::endl;  // not unused
}

}  // namespace mm
