add_library(kdtree
        ../../include/medusa/bits/spatial_search/KDTree.hpp
        ../../include/medusa/bits/spatial_search/KDTree_fwd.hpp
        KDTree.cpp
        ../../include/medusa/bits/spatial_search/KDTreeMutable.hpp
        ../../include/medusa/bits/spatial_search/KDTreeMutable_fwd.hpp
        KDTreeMutable.cpp)

target_link_libraries(kdtree assert types memutils)
register_library(kdtree)


add_library(kdgrid
        ../../include/medusa/bits/spatial_search/Grid_fwd.hpp
        ../../include/medusa/bits/spatial_search/Grid.hpp
        Grid.cpp
        ../../include/medusa/bits/spatial_search/KDGrid_fwd.hpp
        ../../include/medusa/bits/spatial_search/KDGrid.hpp
        KDGrid.cpp)

target_link_libraries(kdgrid assert types memutils)
register_library(kdgrid)
