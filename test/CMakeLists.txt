cmake_minimum_required(VERSION 3.2)
project(seqan3 CXX)

### set googletest
enable_testing()

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++1z -fconcepts")

set(GOOGLETEST_PATH "${CMAKE_SOURCE_DIR}/../vendor/googletest/googletest")
set(SEQAN3_LIBRARIES "pthread" "gtest" "gtest_main")

include(ExternalProject)
ExternalProject_Add(
    googletest
    PREFIX googletest
    GIT_REPOSITORY "https://github.com/google/googletest.git"
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/../vendor/googletest"
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR};-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER};-DBUILD_GTEST=1;-DBUILD_GMOCK=0"
)

add_library (gtest_main STATIC IMPORTED)
add_dependencies(gtest_main googletest)
set_target_properties(gtest_main PROPERTIES IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/lib/libgtest_main.a")

add_library (gtest STATIC IMPORTED)
add_dependencies(gtest gtest_main)
set_target_properties(gtest PROPERTIES IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/lib/libgtest.a")

include_directories("${GOOGLETEST_PATH}/include/")
include_directories("${CMAKE_SOURCE_DIR}/../include/")

add_subdirectory(alphabet)
