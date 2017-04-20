cmake_minimum_required(VERSION 3.2)
project(seqan3_tests CXX)

### set googletest
enable_testing()

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -fconcepts")

# CMAKE_SOURCE_DIR is seqan3/test (the path of this CMakeLists.txt)
set(SEQAN3_ROOT "${CMAKE_SOURCE_DIR}/..")
set(GOOGLETEST_ROOT "${PROJECT_BINARY_DIR}/vendor/googletest")
set(SEQAN3_LIBRARIES "pthread" "gtest" "gtest_main")

find_program(CCACHE_PROGRAM ccache)
# TODO(marehr): Once we require cmake>=3.4 use CMAKE_CXX_COMPILER_LAUNCHER instead
# https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_LAUNCHER.html
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
endif(CCACHE_PROGRAM)

include(ExternalProject)
# TODO(marehr): Once we require cmake>=3.4 add "-DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER};"
# to CMAKE_ARGS
ExternalProject_Add(
    googletest
    PREFIX googletest
    GIT_REPOSITORY "https://github.com/google/googletest.git"
    SOURCE_DIR "${GOOGLETEST_ROOT}"
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR};-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER};-DBUILD_GTEST=1;-DBUILD_GMOCK=0"
    UPDATE_DISCONNECTED yes
)

add_library (gtest_main STATIC IMPORTED)
add_dependencies(gtest_main googletest)
set_target_properties(gtest_main PROPERTIES IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/lib/libgtest_main.a")

add_library (gtest STATIC IMPORTED)
add_dependencies(gtest gtest_main)
set_target_properties(gtest PROPERTIES IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/lib/libgtest.a")

# TODO(marehr): Once we require cmake>=3.7 use BUILDSYSTEM_TARGETS instead
# http://stackoverflow.com/a/36341753 - add gtest as dependency to all tests
macro(add_executable _target)
    _add_executable(${_target} ${ARGN})
    add_dependencies(${_target} gtest)
endmacro()

# define all include folders
include_directories("${GOOGLETEST_ROOT}/googletest/include/")
include_directories("${SEQAN3_ROOT}/include/")
include_directories("${SEQAN3_ROOT}/range-v3/include/")
include_directories("${SEQAN3_ROOT}/sdsl-lite/include/")

# add all test cases
add_subdirectory(alphabet)
add_subdirectory(core)
