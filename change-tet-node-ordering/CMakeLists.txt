
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_STANDARD 11) 
 
PROJECT(change-node-ordering)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(change-node-ordering  main.cpp )
 
target_link_libraries(change-node-ordering ${VTK_LIBRARIES})
