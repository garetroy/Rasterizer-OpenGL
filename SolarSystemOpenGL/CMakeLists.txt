cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(solar)
SET(VTK_DIR /Users/hank/Hartree/VTK/VTK6.0.0)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(solar MACOSX_BUNDLE solar)

if(VTK_LIBRARIES)
  target_link_libraries(solar ${VTK_LIBRARIES})
else()
  target_link_libraries(solar vtkHybrid)
endif()

