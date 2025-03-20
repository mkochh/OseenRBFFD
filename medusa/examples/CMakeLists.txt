include_directories(${HDF5_INCL})
link_directories(${HDF5_LINK})

set_property(GLOBAL PROPERTY ALL_EXAMPLES)
set_property(GLOBAL PROPERTY ALL_RUNNABLE_EXAMPLES)
function (declare_example name)
    set(options RUN)  # if the example is considered short enough to be run during tests
    cmake_parse_arguments(declare_example "${options}" "" "" ${ARGN})

    add_executable(${name} "${name}.cpp")
    target_link_libraries(${name} medusa)
    target_compile_options(${name} PRIVATE "-O3")
    set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    # add to all examples
    get_property(tmp GLOBAL PROPERTY ALL_EXAMPLES)
    list(APPEND tmp ${name})
    set_property(GLOBAL PROPERTY ALL_EXAMPLES ${tmp})
    # add command to run example immediately after compiling from its respective directory
    add_custom_target(${name}_run
            COMMAND echo "Running ${name} ..."
            COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/${name}
            COMMAND echo "Done!"
            DEPENDS ${name}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    if (declare_example_RUN)
        get_property(tmp GLOBAL PROPERTY ALL_RUNNABLE_EXAMPLES)
        list(APPEND tmp "${name}_run")
        set_property(GLOBAL PROPERTY ALL_RUNNABLE_EXAMPLES ${tmp})
    endif()
endfunction(declare_example)

# include all examples in subfolders
add_subdirectory(3d_objects)
add_subdirectory(advection_diffusion)
add_subdirectory(cahnHilliard_equation)
add_subdirectory(coupled_domains)
add_subdirectory(customization)
add_subdirectory(electromagnetic_scattering)
add_subdirectory(ghost_nodes)
add_subdirectory(linear_elasticity)
add_subdirectory(nonNewtonian_fluid)
add_subdirectory(poisson_equation)
add_subdirectory(quantum_mechanics)
add_subdirectory(thermo_fluid)
add_subdirectory(wave_equation)
add_subdirectory(inside_by_oversampling)

# join all examples in a single target
get_property(ALL_EXAMPLES GLOBAL PROPERTY ALL_EXAMPLES)
add_custom_target(examples)
add_dependencies(examples ${ALL_EXAMPLES})

get_property(ALL_RUNNABLE_EXAMPLES GLOBAL PROPERTY ALL_RUNNABLE_EXAMPLES)
add_custom_target(examples_run)
add_dependencies(examples_run ${ALL_RUNNABLE_EXAMPLES})
