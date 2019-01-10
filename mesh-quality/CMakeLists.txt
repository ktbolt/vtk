
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_STANDARD 11) 
 
PROJECT(mesh-quality)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(mesh-quality  main.cxx )
 
target_link_libraries(mesh-quality ${VTK_LIBRARIES})
