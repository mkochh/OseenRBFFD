set(EXECUTABLES 
  oseenHLU
)

foreach(e ${EXECUTABLES})
  add_executable(${e} ${e}.cpp)
  foreach(lib ${LINK_LIBRARIES})
    target_link_libraries(${e} ${lib})    
  endforeach(lib ${LINK_LIBRARIES})

endforeach(e ${EXECUTABLES})
