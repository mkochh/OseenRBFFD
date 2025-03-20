declare_test(end2end FILES
        diffusion_explicit.cpp
        poisson_explicit.cpp
        poisson_implicit.cpp
        complex_numbers.cpp
        coupled_system.cpp
        biharmonic_custom.cpp
        cantilever_beam_implicit.cpp
        poisson_ghost_nodes.cpp
        LINK medusa NO_BASELIB)
target_compile_definitions(end2end_test PUBLIC NDEBUG)