# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project(Hole-Measures)

cmake_minimum_required(VERSION 2.8.10)

# OpenMP
INCLUDE(CheckCXXSourceCompiles)
include_directories (include)
FIND_PACKAGE(OpenMP)
set (CMAKE_REQUIRED_FLAGS ${OpenMP_CXX_FLAGS})
CHECK_CXX_SOURCE_COMPILES("
#include <omp.h>
int main() {
#if (_OPENMP >= 200805 || _MSC_VER >= 1500)
  return 0;
#else
  breaks_on_purpose
#endif
}
" OPENMP_VERSION)

if(OPENMP_VERSION)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
  message("
OpenMP 3.0 not supported by the compiler (${CMAKE_CXX_COMPILER})!
To get optimal performance for the \"chunk\" algorithm, use a C++ compiler with OpenMP 3.0 support (e.g., GCC >=4.4).
To use a different compiler, pass it to cmake in the variable CMAKE_CXX_COMPILER:
 cmake . -DCMAKE_CXX_COMPILER=g++-4.7

")
endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS Qt5 )
if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library. CGAL was not found so the project will not be compiled.")
  return()
endif()


# Boost and its components
find_package( Boost REQUIRED )
if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library. Boost was not found so the project will not be compiled.")
  return()
endif()

# This is necessary for displaying meshes with CGAL::draw()
# See also [1]
if(CGAL_Qt5_FOUND)
  #required to use basic_viewer
  add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
else()
  message(STATUS "This project requires the CGAL_Qt5. CGAL_Qt5 was not found so the project will not be compiled.")
  return()
endif()


# include for local directory
file(GLOB INCLUDES include/*.h)
file(GLOB SOURCES src/*.cpp)
add_library( BaseFiles ${SOURCES} ${INCLUDES} )
# target_precompile_headers( BaseFiles PUBLIC "include/cgal_typedef.h")

# include for local package
include_directories( BEFORE ../include )


# Creating entries for all C++ files with "main" routine
# ##########################################################

# create_single_source_cgal_program( "src/main.cpp" )
# create_single_source_cgal_program( "src/main_resampling.cpp" ) # WIP by Yann-Situ
#
# # [1]
# if(CGAL_Qt5_FOUND)
#   target_link_libraries(main PUBLIC CGAL::CGAL_Qt5)
#   target_link_libraries(main_resampling PUBLIC CGAL::CGAL_Qt5)
# endif()

# Creating entries for all C++ files with "main" routine
# ##########################################################
set(main_link "main_link")
set(main_medial "main_medial")
#set(main_resampling "main_resampling") # WIP
add_executable(main_link "mains/main_link.cpp")
add_executable(main_medial "mains/main_medial.cpp")
#add_executable(main_resampling "mains/main_resampling.cpp") # WIP

# [1]
if(CGAL_Qt5_FOUND)
    target_link_libraries(main_link CGAL::CGAL CGAL::CGAL_Qt5)
    target_link_libraries(main_medial CGAL::CGAL CGAL::CGAL_Qt5)
    #target_link_libraries(main_resampling CGAL::CGAL CGAL::CGAL_Qt5) # WIP

    target_link_libraries(main_link BaseFiles)
    target_link_libraries(main_medial BaseFiles)
    #target_link_libraries(main_resampling BaseFiles) # WIP
endif()
