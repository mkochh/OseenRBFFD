file(GLOB sources *.cpp)
file(GLOB headers *.h)
file(GLOB headers *.hpp)

set(SUBDIRS
  cluster
  sparse
  preconditioning
  solver
  oseen_weights_rbffd
  auxiliaries)

foreach(DIR ${SUBDIRS})
  include_directories(${DIR})
endforeach()

foreach(DIR ${SUBDIRS})
  add_subdirectory(${DIR})
endforeach()

set(extension_sources ${extension_sources} ${sources} PARENT_SCOPE)
set(extension_headers ${extension_sources} ${headers} PARENT_SCOPE)


