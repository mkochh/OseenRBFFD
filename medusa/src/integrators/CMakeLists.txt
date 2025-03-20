add_library(integrators
        ../../include/medusa/bits/integrators/AdamsBashforth_fwd.hpp
        ../../include/medusa/bits/integrators/AdamsBashforth.hpp
        AdamsBashforth.cpp
        ../../include/medusa/bits/integrators/RKExplicit_fwd.hpp
        ../../include/medusa/bits/integrators/RKExplicit.hpp
        RKExplicit.cpp
        )
target_link_libraries(integrators assert types numutils)
register_library(integrators)
