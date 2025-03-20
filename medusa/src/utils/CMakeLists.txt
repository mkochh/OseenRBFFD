add_library(assert
        ../../include/medusa/bits/utils/assert.hpp
        ../../include/medusa/Config.hpp
        assert.cpp)
add_library(memutils
        ../../include/medusa/bits/utils/memutils.hpp
        memutils.cpp)
add_library(numutils
        ../../include/medusa/bits/utils/numutils.hpp
        numutils.cpp)
target_link_libraries(numutils assert)

add_library(randutils
        ../../include/medusa/bits/utils/randutils.hpp
        randutils.cpp)
target_link_libraries(randutils assert)

add_library(stdtypesutils
        ../../include/medusa/bits/utils/stdtypesutils.hpp
        stdtypesutils.cpp)
target_link_libraries(stdtypesutils assert)

add_library(timer
        ../../include/medusa/bits/utils/Timer.hpp
        Timer.cpp)
target_link_libraries(timer assert)

add_library(stopwatch
        ../../include/medusa/bits/utils/Stopwatch.hpp
        Stopwatch.cpp)
target_link_libraries(stopwatch assert)

register_library(assert)
register_library(memutils)
register_library(stdtypesutils)
register_library(randutils)
register_library(timer)
register_library(stopwatch)
register_library(numutils)