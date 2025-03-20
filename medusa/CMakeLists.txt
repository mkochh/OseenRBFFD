cmake_minimum_required(VERSION 3.25)

# Setting compilers (must be before first project command)
# this has to be set in the main CMakeLists file of the project

# SET(CMAKE_C_COMPILER /usr/bin/gcc)
# SET(CMAKE_CXX_COMPILER /usr/bin/g++)

# SET(CMAKE_C_COMPILER /usr/bin/clang)
# SET(CMAKE_CXX_COMPILER /usr/bin/clang++)

# SET(CMAKE_C_COMPILER /usr/local/opt/llvm/bin/clang)
# SET(CMAKE_CXX_COMPILER /usr/local/opt/llvm/bin/clang++)

# SET(CMAKE_C_COMPILER /opt/intel/bin/icc)
# SET(CMAKE_CXX_COMPILER /opt/intel/bin/icpc)

project(medusa)
set(MEDUSA_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
include(${MEDUSA_ROOT}/scripts/medusa.cmake)

# Do not include unnecessary compile flags if added as a subdirectory.
if(${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wfloat-conversion -Wno-deprecated-copy -Wno-maybe-uninitialized -pedantic -O3 -std=c++14 ${CMAKE_CXX_FLAGS}")
else()
    set(CMAKE_CXX_FLAGS "-std=c++17 ${CMAKE_CXX_FLAGS}")
endif()

include_directories(SYSTEM ${MEDUSA_ROOT}/include/)
add_subdirectory(${MEDUSA_ROOT}/src/)

# Do not include unnecessary targets if added as a subdirectory.
if(${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    add_subdirectory(${MEDUSA_ROOT}/test/)
    add_subdirectory(${MEDUSA_ROOT}/examples/)
    add_subdirectory(${MEDUSA_ROOT}/docs/)
endif()
