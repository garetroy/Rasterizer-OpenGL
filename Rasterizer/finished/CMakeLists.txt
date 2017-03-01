cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(main)
SET(VTK_DIR /Users/hank/Hartree/VTK/VTK6.0.0)

file(GLOB SRCFILES
    "*.h"
    "*.cxx"
)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(main MACOSX_BUNDLE ${SRCFILES} main)

if(VTK_LIBRARIES)
  target_link_libraries(main ${VTK_LIBRARIES})
else()
  target_link_libraries(main vtkHybrid)
endif()

