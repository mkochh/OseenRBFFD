include_directories(${HDF5_INCL})
link_directories(${HDF5_LINK})

set_property(GLOBAL PROPERTY ALL_SOURCES)
set_property(GLOBAL PROPERTY ALL_LIBS)
function(register_library target)
    # add sources to list of all sources
    get_target_property(src ${target} SOURCES)
    get_property(tmp GLOBAL PROPERTY ALL_SOURCES)
    foreach(arg ${src})
        list(APPEND tmp "${CMAKE_CURRENT_SOURCE_DIR}/${arg}")
    endforeach()
    set_property(GLOBAL PROPERTY ALL_SOURCES ${tmp})
    # add to list of all libs
    get_property(tmp GLOBAL PROPERTY ALL_LIBS)
    list(APPEND tmp ${target})
    set_property(GLOBAL PROPERTY ALL_LIBS ${tmp})
endfunction()

add_subdirectory(approximations/)
add_subdirectory(interpolants/)
add_subdirectory(io/)
add_subdirectory(domains/)
add_subdirectory(types/)
add_subdirectory(utils/)
add_subdirectory(spatial_search/)
add_subdirectory(operators/)
add_subdirectory(integrators/)

# Library representing all libraries together, as a shortcut for linking.
# Note that the produced .a file does not contain all functions -- use medusa_standalone for that.
add_library(medusa ${MEDUSA_ROOT}/src/types/Vec.cpp)
get_property(ALL_LIBS GLOBAL PROPERTY ALL_LIBS)
target_link_libraries(medusa ${ALL_LIBS})

# Build all aources into one large library that library can be shipped independently.
# The resulting medusa_standalone.a file can be shipped independently only with files in include/.
get_property(ALL_SOURCES GLOBAL PROPERTY ALL_SOURCES)
add_library(medusa_standalone ${ALL_SOURCES})
