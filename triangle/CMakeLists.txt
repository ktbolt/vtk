
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_STANDARD 11) 
 
PROJECT(tri)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(tri main.cpp )
 
target_link_libraries(tri ${VTK_LIBRARIES})
