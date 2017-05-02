cmake_minimum_required(VERSION 3.2)
project(seqan3_tests CXX)

### set googletest
enable_testing()

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -fconcepts")

# CMAKE_SOURCE_DIR is seqan3/test (the path of this CMakeLists.txt)
set(SEQAN3_ROOT "${CMAKE_SOURCE_DIR}/..")
set(GOOGLETEST_ROOT "${PROJECT_BINARY_DIR}/vendor/googletest")

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

macro(seqan3_test target_cpp)
    # $target_cpp = "pod_tuple_test.cpp"
    #   * will register the global TARGET name "pod_tuple_test" and
    #   * will register the test case name "core_pod_tuple" if
    #     pod_tuple_test.cpp is in test/core/
    #
    # NOTE(marehr): ".+/test/" REGEX is greedy, that means
    # /test/test/test/hello_test.cpp will result in an empty `test_path`
    string(REGEX REPLACE "_test.cpp$" "" target_name ${target_cpp})
    string(REGEX REPLACE ".+/test/" "" test_path ${CMAKE_CURRENT_LIST_DIR})
    string(REGEX REPLACE "/" "_" test_path ${test_path})

    # add_executable(pod_tuple_test pod_tuple_test.cpp)
    # target_link_libraries(pod_tuple_test "pthread" "gtest" "gtest_main")
    # add_test(NAME "core_pod_tuple" COMMAND ${target_name}_test)
    # add_dependencies(pod_tuple_test gtest)
    add_executable(${target_name}_test ${target_cpp})
    target_link_libraries(${target_name}_test "pthread" "gtest" "gtest_main")
    add_test(NAME "${test_path}_${target_name}" COMMAND ${target_name}_test)
    add_dependencies(${target_name}_test gtest)
endmacro()

# define all include folders
include_directories("${GOOGLETEST_ROOT}/googletest/include/")
include_directories("${SEQAN3_ROOT}/include/")
include_directories("${SEQAN3_ROOT}/range-v3/include/")
include_directories("${SEQAN3_ROOT}/sdsl-lite/include/")

macro (add_subdirectories)
    file (GLOB ENTRIES
          RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/[!.]*)

    foreach (ENTRY ${ENTRIES})
        if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY})
            if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY}/CMakeLists.txt)
                add_subdirectory(${ENTRY})
            endif ()
        endif ()
    endforeach ()
endmacro ()

add_subdirectories()
