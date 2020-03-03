cmake_minimum_required(VERSION 2.6)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

#set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/.. )

project (high-order-lz77)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

# hopscotch_map
ExternalProject_Add(hopscotch_map
  GIT_REPOSITORY "https://github.com/Tessil/hopscotch-map.git"
  BUILD_IN_SOURCE TRUE
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${INSTALL_DIR} # TODO ADD static build flag
  UPDATE_COMMAND ""
  INSTALL_COMMAND ""
  BUILD_COMMAND ""
  CONFIGURE_COMMAND "")
ExternalProject_Get_property(hopscotch_map INSTALL_DIR)
set(hopscotch_map_INCLUDE "${INSTALL_DIR}/src/hopscotch_map/include/")

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/DYNAMIC/include)
include_directories(${INSTALL_DIR}/src/hopscotch_map/include)

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS "--std=c++11")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g")
set(CMAKE_CXX_FLAGS_RELEASE "-g -ggdb -Ofast -fstrict-aliasing -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")

add_executable(ho-lz77 ho-lz77.cpp)
add_dependencies(ho-lz77 hopscotch_map)



