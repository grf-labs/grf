project(ranger_test)
cmake_minimum_required(VERSION 2.0)

## ======================================================================================##
## Compiler flags
## ======================================================================================##
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
set_target_properties(${TARGET} PROPERTIES COMPILE_FLAGS -pthread)

## ======================================================================================##
## Dependencies
## ======================================================================================##
include_directories(../source/src ../source/src/utility ../source/src/Forest ../source/src/Tree)
add_subdirectory(../source ranger_build)

## ======================================================================================##
## Subdirectories and source files
## ======================================================================================##
file(GLOB TEST_SOURCES *.cpp)
file(GLOB_RECURSE RG_SOURCES ../source/src/*.cpp)

## Remove main file from ranger
get_filename_component(MAIN_FILE . ABSOLUTE)
SET(MAIN_FILE "${MAIN_FILE}/../source/src/main.cpp")
list(REMOVE_ITEM RG_SOURCES "${MAIN_FILE}")
##

set(SOURCES "${RG_SOURCES};${TEST_SOURCES}")

## ======================================================================================##
## Google test
## ======================================================================================##
ADD_SUBDIRECTORY (gtest-1.7.0)
enable_testing()
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
add_test(runUnitTests runUnitTests)

## ======================================================================================##
## Executable
## ======================================================================================##
add_executable(runUnitTests ${SOURCES})
add_dependencies(runUnitTests ranger)
target_link_libraries(runUnitTests gtest gtest_main)

