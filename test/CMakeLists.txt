# ============================================================================
#                  SeqAn - The Library for Sequence Analysis
# ============================================================================
#
# Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
# Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Knut Reinert or the FU Berlin nor the names of
#       its contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
# ============================================================================

cmake_minimum_required (VERSION 3.2)

project (seqan3_tests CXX)

include (ExternalProject)

# ----------------------------------------------------------------------------
# Define paths
# ----------------------------------------------------------------------------

# CMAKE_SOURCE_DIR is seqan3/test (the path of this CMakeLists.txt)
set(SEQAN3_ROOT "${CMAKE_SOURCE_DIR}/..")

# SeqAn include
set (CMAKE_INCLUDE_PATH ${SEQAN3_ROOT}/include ${CMAKE_INCLUDE_PATH})

# SeqAn cmake module
set (CMAKE_PREFIX_PATH ${SEQAN3_ROOT}/build_system ${CMAKE_PREFIX_PATH})
set (CMAKE_MODULE_PATH ${SEQAN3_ROOT}/build_system ${CMAKE_MODULE_PATH})

# Auto-add dependencies if present in same dir (otherwise will be searched globally)
foreach (_PATH
         "${SEQAN3_ROOT}/range-v3/include/"
         "${SEQAN3_ROOT}/sdsl-lite/include/"
         "${SEQAN3_ROOT}/cereal/include/"
         "${SEQAN3_ROOT}/lemon/")
    if (EXISTS ${_PATH})
        set (CMAKE_INCLUDE_PATH ${_PATH} ${CMAKE_INCLUDE_PATH})
    endif ()
endforeach ()

include_directories (${CMAKE_INCLUDE_PATH})

# ----------------------------------------------------------------------------
# Load SeqAn module
# ----------------------------------------------------------------------------

find_package (SeqAn3 REQUIRED)

# add include dir of test directory for some global test utility functions.
set (SEQAN3_INCLUDE_DIRS "${SEQAN3_INCLUDE_DIRS};${CMAKE_SOURCE_DIR}/include")

# import package's information (libraries are added with targets below)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN3_CXX_FLAGS}")
include_directories (${SEQAN3_INCLUDE_DIRS})
add_definitions (${SEQAN3_DEFINITIONS})

# add warning levels
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Werror -Wall -Wextra")

# ----------------------------------------------------------------------------
# Use CCache if possible
# ----------------------------------------------------------------------------

find_program (CCACHE_PROGRAM ccache)
# TODO(marehr): Once we require cmake>=3.4 use CMAKE_CXX_COMPILER_LAUNCHER instead
# https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_LAUNCHER.html
if (CCACHE_PROGRAM)
    set_property (GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
    set_property (GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
endif (CCACHE_PROGRAM)

# ----------------------------------------------------------------------------
# Add GTest
# ----------------------------------------------------------------------------

enable_testing ()

set (GOOGLETEST_ROOT "${PROJECT_BINARY_DIR}/vendor/googletest")
include_directories ("${GOOGLETEST_ROOT}/googletest/include/")

# TODO(marehr): Once we require cmake>=3.4 add "-DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER};"
# to CMAKE_ARGS
ExternalProject_Add(
    googletest
    PREFIX googletest
    GIT_REPOSITORY "https://github.com/google/googletest.git"
    GIT_TAG "15392f1a38fa0b8c3f13a9732e94b209069efa1c"
    SOURCE_DIR "${GOOGLETEST_ROOT}"
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR};-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER};-DBUILD_GTEST=1;-DBUILD_GMOCK=0"
    UPDATE_DISCONNECTED yes
)

add_library (gtest_main STATIC IMPORTED)
add_dependencies (gtest_main googletest)
set_target_properties (gtest_main PROPERTIES IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/lib/libgtest_main.a")

add_library (gtest STATIC IMPORTED)
add_dependencies (gtest gtest_main)
set_target_properties (gtest PROPERTIES IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/lib/libgtest.a")

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

    add_executable(${target_name}_test ${target_cpp})
    target_link_libraries(${target_name}_test ${SEQAN3_LIBRARIES} "gtest" "gtest_main")
    add_test(NAME "${test_path}_${target_name}" COMMAND ${target_name}_test)
    add_dependencies(${target_name}_test gtest)
endmacro()

# ----------------------------------------------------------------------------
# Add the tests
# ----------------------------------------------------------------------------

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
