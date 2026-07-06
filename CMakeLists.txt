
# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.31)

project(Hole-Measures)

if (${APPLE})
    message("Apple")
    set(CMAKE_PREFIX_PATH "/opt/homebrew/opt")
    set (CMAKE_CXX_STANDARD 17)
    set(CGAL_DIR "/opt/homebrew/opt/cgal/lib/cmake/CGAL")
    include_directories(/opt/homebrew/opt/cgal/include)
    include_directories(/opt/homebrew/opt/boost/include)
endif()

INCLUDE(CheckCXXSourceCompiles)

include_directories (include)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS)
find_package( Boost REQUIRED )

# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
  create_single_source_cgal_program("${cppfile}")
endforeach()
