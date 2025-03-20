add_library(hdf
        ../../include/medusa/bits/io/HDF.hpp
        ../../include/medusa/bits/io/HDF_fwd.hpp
        ../../include/medusa/bits/io/HDF_Eigen.hpp
        HDF.cpp)
target_link_libraries(hdf assert hdf5)

add_library(xml
        ../../include/medusa/bits/io/XML.hpp
        ../../include/medusa/bits/io/XML_fwd.hpp
        XML.cpp)
target_link_libraries(xml assert)

add_library(csv
        ../../include/medusa/bits/io/CSV.hpp
        ../../include/medusa/bits/io/CSV_fwd.hpp
        ../../include/medusa/bits/io/CSV_Eigen.hpp
        CSV.cpp)
target_link_libraries(csv assert)

add_library(stl
        ../../include/medusa/bits/io/STL.hpp
        ../../include/medusa/bits/io/STL_fwd.hpp
        STL.cpp)
target_link_libraries(stl assert)

add_library(ioformats ../../include/medusa/bits/io/ioformat.hpp ioformats.cpp)

register_library(hdf)
register_library(csv)
register_library(xml)
register_library(stl)
register_library(ioformats)
