cmake_minimum_required(VERSION 3.25)

project("MedusaH2LibExtended")

# Set build type to release, if not specified.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Add math library to linked libraries
set(LINK_LIBRARIES m)

# Add definition for H2Lib quick exit at rank 0.
add_definitions(-DHARITH_RKMATRIX_QUICK_EXIT)

# Add BLAS library to linked libraries
set(BLA_VENDOR OpenBLAS)
find_package(BLAS)
if(BLAS_FOUND)
  message("-- " ${BLAS_LIBRARIES})
  add_definitions(-DUSE_BLAS)
  set(LINK_LIBRARIES ${LINK_LIBRARIES} ${BLAS_LIBRARIES})
else()
  message("-- BLAS not found")
endif(BLAS_FOUND)

find_package(OpenMP)
if (OPENMP_FOUND)
  message("-- OpenMP found")
  add_definitions(-DUSE_OPENMP)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
  message("-- OpenMP not found")
endif()

# Build extended H2Lib
add_subdirectory(H2Lib)
include_directories(H2Lib)
link_directories(H2Lib)

# Add H2Lib to linked libraries
set(LINK_LIBRARIES ${LINK_LIBRARIES} h2lib)

add_subdirectory(${CMAKE_SOURCE_DIR}/medusa medusa)
include_directories(${CMAKE_SOURCE_DIR}/medusa/include)
set(LINK_LIBRARIES ${LINK_LIBRARIES} medusa)

set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
find_package(CGAL REQUIRED)
if(CGAL_FOUND)
  message("-- Using CGAL Version ${CGAL_VERSION}")
  if(CGAL_VERSION_MAJOR LESS 6)
    add_definitions(-DUSE_CGAL_LESS_THAN_6)
  else()
    add_definitions(-DUSE_CGAL6)
  endif()
endif()

# Build Extensions
# Set flags for debugging or optimization.
set(CMAKE_CXX_FLAGS "-std=c++17 -fopenmp -Wall -Wextra -Wno-unused-function -Wno-maybe-uninitialized -Wno-deprecated-copy -Wno-deprecated-declarations -Wno-overloaded-virtual -g -O3")
set(CMAKE_C_FLAGS "-fopenmp -Wall -Wextra -Wno-unused-function -Wno-maybe-uninitialized -Wno-deprecated-copy -Wno-deprecated-declarations -Wno-overloaded-virtual -g -O3")

set(extension_headers)
set(extension_sources)
add_subdirectory(Extensions)
include_directories(Extensions)

add_library(h2libext ${extension_sources})

foreach(lib ${LINK_LIBRARIES})
  target_link_libraries(h2libext ${lib})
endforeach(lib ${LINK_LIBRARIES})

set(LINK_LIBRARIES ${LINK_LIBRARIES} h2libext)

# Build Tests
add_subdirectory(Tests)
