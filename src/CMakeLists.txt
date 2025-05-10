list(APPEND polyhedron_headers ${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp)
list(APPEND polyhedron_sources ${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp)

list(APPEND polyhedron_headers ${CMAKE_CURRENT_SOURCE_DIR}/Polyhedron.hpp)

list(APPEND polyhedron_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(polyhedron_sources ${polyhedron_sources} PARENT_SCOPE)
set(polyhedron_headers ${polyhedron_headers} PARENT_SCOPE)
set(polyhedron_includes ${polyhedron_includes} PARENT_SCOPE)