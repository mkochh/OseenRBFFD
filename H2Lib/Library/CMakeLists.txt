# Add all source files to library files
set(H2LIB_SRCS
    aca.c
    amatrix.c
    avector.c
    basic.c
    bem2d.c
    bem3d.c
    block.c
    clusterbasis.c
    clustergeometry.c
    cluster.c
    clusteroperator.c
    curve2d.c
    dblock.c
    dclusterbasis.c
    dcluster.c
    dclusteroperator.c
    ddcluster.c
    dh2compression.c
    dh2matrix.c
    duniform.c
    eigensolvers.c
    factorizations.c
    gaussquad.c
    h2arith.c
    h2compression.c
    h2matrix.c
    h2update.c
    harith2.c
    harith.c
    hcoarsen.c
    helmholtzbem3d.c
    helmholtzoclbem3d.c
    hmatrix.c
    krylov.c
    krylovsolvers.c
    laplacebem2d.c
    laplacebem3d.c
    laplaceoclbem3d.c
    macrosurface3d.c
    matrixnorms.c
    oclbem3d.c
    opencl.c
    parameters.c
    realavector.c
    rkmatrix.c
    settings.c
    singquad1d.c
    singquad2d.c
    sparsematrix.c
    sparsepattern.c
    surface3d.c
    tet3d.c
    tet3dp1.c
    tet3drt0.c
    tri2d.c
    tri2dp1.c
    tri2drt0.c
    truncation.c
    uniform.c
    visualize.c)

set(H2LIB_HDRS
    aca.h
    amatrix.h
    avector.h
    basic.h
    bem2d.h
    bem3d.h
    blas.h
    block.h
    clsettings.h
    clusterbasis.h
    clustergeometry.h
    cluster.h
    clusteroperator.h
    curve2d.h
    dblock.h
    dclusterbasis.h
    dcluster.h
    dclusteroperator.h
    ddcluster.h
    dh2compression.h
    dh2matrix.h
    duniform.h
    eigensolvers.h
    factorizations.h
    gaussquad.h
    h2arith.h
    h2compression.h
    h2matrix.h
    h2update.h
    harith2.h
    harith.h
    hcoarsen.h
    helmholtzbem3d.h
    helmholtzoclbem3d.h
    hmatrix.h
    krylov.h
    krylovsolvers.h
    laplacebem2d.h
    laplacebem3d.h
    laplaceoclbem3d.h
    macrosurface3d.h
    matrixnorms.h
    oclbem3d.h
    opencl.h
    parameters.h
    realavector.h
    rkmatrix.h
    settings.h
    simd_avx.h
    simd.h
    simd_sse2.h
    singquad1d.h
    singquad2d.h
    sparsematrix.h
    sparsepattern.h
    surface3d.h
    tet3d.h
    tet3dp1.h
    tet3drt0.h
    tri2d.h
    tri2dp1.h
    tri2drt0.h
    truncation.h
    uniform.h
    visualize.h)

if(USE_CAIRO)
    include_directories(${CAIRO_INCLUDE_DIR})
endif()

set(sources)

foreach(src ${H2LIB_SRCS})
    set(sources ${sources} ${CMAKE_CURRENT_SOURCE_DIR}/${src})
endforeach(src ${H2LIB_SRC})

set(headers)

foreach(hdr ${H2LIB_HDRS})
    set(headers ${headers} ${CMAKE_CURRENT_SOURCE_DIR}/${hdr})
endforeach(hdr ${H2LIB_HDR})

set(library_sources ${library_sources} ${sources} PARENT_SCOPE)
set(library_headers ${library_headers} ${headers} PARENT_SCOPE)