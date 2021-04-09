
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(ImageReslice)

find_package(VTK COMPONENTS 
  vtkCommonCore
  vtkCommonDataModel
  vtkImagingCore
  vtkInteractionStyle
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingOpenGL2 QUIET)
if (NOT VTK_FOUND)
  message("Skipping ImageReslice: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(ImageReslice ImageReslice.cxx )
  target_link_libraries(ImageReslice PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  add_executable(ImageReslice ImageReslice.cxx )
  #add_executable(ImageReslice MACOSX_BUNDLE ImageReslice.cxx )
  target_link_libraries(ImageReslice PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS ImageReslice
    MODULES ${VTK_LIBRARIES}
    )
endif () 
