
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(ImageSlice)

find_package(VTK COMPONENTS 
  vtkCommonCore
  vtkCommonDataModel
  vtkInteractionStyle
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingImage
  vtkRenderingOpenGL2 QUIET)
if (NOT VTK_FOUND)
  message("Skipping ImageSlice: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(ImageSlice ImageSlice.cxx )
  target_link_libraries(ImageSlice PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  add_executable(ImageSlice ImageSlice.cxx )
  target_link_libraries(ImageSlice PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS ImageSlice
    MODULES ${VTK_LIBRARIES}
    )
endif () 
