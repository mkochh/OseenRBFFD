declare_test(domains FILES
        discretization_helpers_test.cpp
        discretization_helpers_advanced_test.cpp
        DomainDiscretization_test.cpp
        DomainShape_test.cpp
        BoxShape_test.cpp
        BallShape_test.cpp
        PolygonShape_test.cpp
        PolytopeShape_test.cpp
        RotatedShape_test.cpp
        UnknownShape_test.cpp
        ShapeUnion_test.cpp
        ShapeDifference_test.cpp
        STLShape_test.cpp
        NURBSShape_test.cpp
        TranslatedShape_test.cpp
        LINK hdf stl)

declare_test(domain_engines FILES
        GeneralFill_test.cpp
        GeneralSurfaceFill_test.cpp
        compute_normal_test.cpp
        cad_helpers_test.cpp
        NURBSPatch_test.cpp
        GrainDropFill_test.cpp
        GridFill_test.cpp
        FindClosest_test.cpp
        FindBalancedSupport_test.cpp
        HalfLinksRefine_test.cpp
        BasicRelax_test.cpp)

if (CGAL_FOUND)
    declare_test(domains_cgal FILES
            PolyhedronShape_test.cpp)
endif()