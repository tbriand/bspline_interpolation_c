cmake_minimum_required(VERSION 2.6)

# FindGSL.cmake not present before CMake 3.2, use provided one if necessary
if(CMAKE_VERSION VERSION_LESS "3.2")
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake")
endif()

set(GSL_FIND_QUIETLY TRUE)
find_package(GSL)
if(GSL_FOUND)
  include_directories("${GSL_INCLUDE_DIRS}")
  add_definitions(-DGSL_SUPPORT)
endif()

# IIO
if(NOT TARGET IIOLIB)
   include_directories("iio")
   add_subdirectory(iio)
endif()

# Enable C99
if(CMAKE_VERSION VERSION_LESS "3.1")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
else()
  set(CMAKE_C_STANDARD 99)
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")

add_library(Splinter bspline.c bspline.h splinter.c splinter.h)
target_link_libraries(Splinter ${GSL_LIBRARIES} m)

add_executable(bspline bspline_main.c splinter_transform.c homography_tools.c)
target_link_libraries(bspline IIOLIB Splinter)

add_executable (hom4p hom4p.c homography_tools.c)
target_link_libraries(hom4p m)

if(GSL_FOUND)
  add_executable (compute_bspline compute_bspline.c bspline.c)
  target_link_libraries(compute_bspline ${GSL_LIBRARIES} m)
endif()
