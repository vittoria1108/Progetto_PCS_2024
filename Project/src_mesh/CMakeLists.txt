list(APPEND src_mesh_headers ${CMAKE_CURRENT_SOURCE_DIR}/PolygonalMesh.hpp)
list(APPEND src_mesh_sources ${CMAKE_CURRENT_SOURCE_DIR}/PolygonalMesh.cpp)
list(APPEND src_mesh_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(src_mesh_sources ${src_mesh_sources} PARENT_SCOPE)
set(src_mesh_headers ${src_mesh_headers} PARENT_SCOPE)
set(src_mesh_includes ${src_mesh_includes} PARENT_SCOPE)
