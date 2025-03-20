include_directories(${HDF5_INCL})
link_directories(${HDF5_LINK})

set(GTEST_ROOT gtest-1.7.0)
set(GTEST_LIBS gtest gtest_main)
include_directories(SYSTEM ${GTEST_ROOT}/include)
add_subdirectory(${GTEST_ROOT})

# Usage:
#
# declate_test(component_name
#     FILES test_file1.cpp test_file2.cpp
#     [LINK additional_lib1 aditional_lib2]
#     [NO_BASELIB])
#
# Test is automatically linked to google test libs and its component lib.
# Addidional libraries needed only by the test can be added after LINK keyword.
# NO_BASELIB option can be specifies for tests that do not test a library with the same name.
#
# Using environment variable GF, filters can be passed to Google test.
set_property(GLOBAL PROPERTY ALL_TESTS)
function(declare_test test_name)
    # parse args
    set(options NO_BASELIB EXCLUDE_FROM_ALL)
    set(multiValueArgs FILES LINK)
    cmake_parse_arguments(declare_test "${options}" "" "${multiValueArgs}" ${ARGN})
    add_executable(${test_name}_test ${declare_test_FILES})
    # link test to its appropriate library
    if (NOT declare_test_NO_BASELIB)
        target_link_libraries(${test_name}_test ${test_name})
    endif()
    # link test to other libraries
    target_link_libraries(${test_name}_test ${GTEST_LIBS} ${declare_test_LINK})
    # add command to run test immediately after compiling from MEDUSA_ROOT directory
    # with enabled support for filters
    add_custom_target(${test_name}_run_tests
            COMMAND ${MEDUSA_ROOT}/bin/${test_name}_test --gtest_filter=*$(GF)*
            WORKING_DIRECTORY ${MEDUSA_ROOT})
    add_dependencies(${test_name}_run_tests ${test_name}_test)
    # add this test to list of all tests (unless specified otherwise)
    if (NOT declare_test_EXCLUDE_FROM_ALL)
        get_property(tmp GLOBAL PROPERTY ALL_TESTS)
        foreach(arg ${declare_test_FILES})
            list(APPEND tmp "${CMAKE_CURRENT_SOURCE_DIR}/${arg}")
        endforeach()
        set_property(GLOBAL PROPERTY ALL_TESTS ${tmp})
    endif()
endfunction(declare_test)

add_subdirectory(types/)
add_subdirectory(approximations/)
add_subdirectory(interpolants/)
add_subdirectory(utils/)
add_subdirectory(spatial_search/)
add_subdirectory(domains/)
add_subdirectory(operators/)
add_subdirectory(io/)
add_subdirectory(integrators/)
add_subdirectory(end2end/)

# all tests together
get_property(ALL_TESTS GLOBAL PROPERTY ALL_TESTS)
add_executable(medusa_tests ${ALL_TESTS})
target_link_libraries(medusa_tests medusa ${GTEST_LIBS})
add_custom_target(medusa_run_tests DEPENDS medusa_tests
        COMMAND ${MEDUSA_ROOT}/bin/medusa_tests
        WORKING_DIRECTORY ${MEDUSA_ROOT})
