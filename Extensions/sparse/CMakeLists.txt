file(GLOB sources *.cpp)
file(GLOB headers *.h)

set(extension_sources  ${extension_sources} ${sources} PARENT_SCOPE)
set(extension_headers  ${extension_sources} ${headers} PARENT_SCOPE)