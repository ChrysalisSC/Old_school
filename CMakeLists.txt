cmake_minimum_required(VERSION 3.10)

PROJECT(astro)
SET(VTK_DIR /home/colton/Desktop/Test/build)



find_package(VTK REQUIRED)

add_executable(astro astro.cxx)

#SET(CMAKE_CXX_LINK_FLAGS "-framework OpenGL")
if(VTK_LIBRARIES)
  target_link_libraries(astro ${VTK_LIBRARIES})
  set_target_properties(astro PROPERTIES LINKER_LANGUAGE CXX)
else()
  target_link_libraries(astro vtkHybrid)
endif()

vtk_module_autoinit(
  TARGETS astro
  MODULES ${VTK_LIBRARIES})
