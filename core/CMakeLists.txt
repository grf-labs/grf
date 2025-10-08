
project(grf)
cmake_minimum_required(VERSION 3.5)
## ======================================================================================##
## Compiler flags
## ======================================================================================##
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2")
if(MSVC)
  # Decrease warning level for verbose MSVC
  add_compile_options(/W3)
endif()
## ======================================================================================##
## Linker flags
## ======================================================================================##
if(NOT MSVC)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")
endif()
## ======================================================================================##
## Subdirectories and source files
## ======================================================================================##
include_directories(src test third_party)
file(GLOB_RECURSE SOURCES src/*.cpp test/*.cpp)
## ======================================================================================##
## Debug and release targets
## ======================================================================================##
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )
ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )
## ======================================================================================##
## Executable
## ======================================================================================##
add_executable(grf ${SOURCES})
