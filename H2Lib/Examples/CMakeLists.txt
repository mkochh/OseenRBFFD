set(EXAMPLES
  example_amatrix_bem3d
  example_dh2matrix
  example_fem_p1_harith
  example_fem_rt0_harith
  example_h2matrix_bem3d
  example_hmatrix_bem3d
  example_visualize)

foreach(e ${EXAMPLES})
  add_executable(${e} ${e}.c)
  target_link_libraries(${e} h2lib m)

  if(USE_BLAS)
    target_link_libraries(${e} blas lapack)
  endif()

  if(USE_CAIRO)
    target_link_libraries(${e} cairo)
  endif()
endforeach()
