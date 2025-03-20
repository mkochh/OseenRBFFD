# Build library
# set(CMAKE_C_FLAGS "-Wall")

include_directories(Library)

set(library_sources)
set(library_headers)

add_subdirectory(Library)

add_library(h2lib ${library_sources})

foreach(lib ${LINK_LIBRARIES})
  target_link_libraries(h2lib ${lib})
endforeach(lib ${LINK_LIBRARIES})

# Add tests from Tests directory
add_subdirectory(Tests)

# Add Examples
add_subdirectory(Examples)