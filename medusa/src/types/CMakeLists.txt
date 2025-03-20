add_library(types
        ../../include/medusa/bits/types/traits.hpp
        ../../include/medusa/bits/types/ScalarField_fwd.hpp
        ../../include/medusa/bits/types/ScalarField.hpp
        ScalarField.cpp
        ../../include/medusa/bits/types/Range_fwd.hpp
        ../../include/medusa/bits/types/Range.hpp
        Range.cpp
        ../../include/medusa/bits/types/Vec_fwd.hpp
        ../../include/medusa/bits/types/Vec.hpp
        Vec.cpp
        ../../include/medusa/bits/types/VectorField_fwd.hpp
        ../../include/medusa/bits/types/VectorField.hpp
        VectorField.cpp)
target_link_libraries(types assert)

register_library(types)
