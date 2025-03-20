set(HDF5_INCL "/usr/include/hdf5/serial/" "/usr/include/x86_64-linux-gnu/hdf5/serial/" "/usr/local/opt/hdf5/include/")
set(HDF5_LINK "/usr/lib/x86_64-linux-gnu/hdf5/serial" "/usr/local/opt/hdf5/lib/")
try_compile(HDF5_TRY_COMPILE
    ${CMAKE_BINARY_DIR}
    ${MEDUSA_ROOT}/test/test_hdf_compile.cpp
    LINK_LIBRARIES hdf5
    CMAKE_FLAGS
      "-DINCLUDE_DIRECTORIES=${HDF5_INCL}"
      "-DLINK_DIRECTORIES=${HDF5_LINK}"
    OUTPUT_VARIABLE HDF5_TEST_OUTPUT)

if (NOT HDF5_TRY_COMPILE)
    message("Compile log: ${HDF5_TEST_OUTPUT}")
    message(FATAL_ERROR "HDF5 sample failed to compile. See errors above.")
else()
    message("-- HDF5 sample compiled sucessfully.")
endif()