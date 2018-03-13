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

project(seqan3_performance_test CXX)

# ----------------------------------------------------------------------------
# Define paths
# ----------------------------------------------------------------------------

# CMAKE_SOURCE_DIR is seqan3/test (the path of this CMakeLists.txt)
set (SEQAN3_ROOT "${CMAKE_SOURCE_DIR}/../..")

# SeqAn cmake module
set (CMAKE_PREFIX_PATH ${SEQAN3_ROOT}/build_system ${CMAKE_PREFIX_PATH})
set (CMAKE_MODULE_PATH ${SEQAN3_ROOT}/build_system ${CMAKE_MODULE_PATH})

# ----------------------------------------------------------------------------
# Load SeqAn module
# ----------------------------------------------------------------------------

find_package (SeqAn3 REQUIRED)

# add include dir of test directory for some global test utility functions.
set (SEQAN3_INCLUDE_DIRS ${SEQAN3_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/../include)

# import package's information (libraries are added with targets below)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN3_CXX_FLAGS}")
include_directories (${SEQAN3_INCLUDE_DIRS})
add_definitions (${SEQAN3_DEFINITIONS})

# add warning levels
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Werror -Wall -Wextra")

# CMAKE_SOURCE_DIR is seqan3/benchmark (the path of this CMakeLists.txt)
set(SEQAN3_ROOT "${CMAKE_SOURCE_DIR}/..")
set(BENCHMARK_ROOT "${PROJECT_BINARY_DIR}/vendor/benchmark")

# ----------------------------------------------------------------------------
# Use CCache if possible
# ----------------------------------------------------------------------------

find_program(CCACHE_PROGRAM ccache)
# TODO(marehr): Once we require cmake>=3.4 use CMAKE_CXX_COMPILER_LAUNCHER instead
# https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_LAUNCHER.html
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
endif(CCACHE_PROGRAM)

# ----------------------------------------------------------------------------
# Add google/benchmark
# ----------------------------------------------------------------------------

include(ExternalProject)
# TODO(marehr): Once we require cmake>=3.4 add "-DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER};"
# to CMAKE_ARGS
set(
    google_benchmark_args
    "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}"
    "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
    "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
    "-DBENCHMARK_ENABLE_TESTING=false"
    # "-DBENCHMARK_ENABLE_LTO=true"
)
ExternalProject_Add(
    google_benchmark
    PREFIX google_benchmark
    GIT_REPOSITORY "https://github.com/google/benchmark.git"
    SOURCE_DIR "${BENCHMARK_ROOT}"
    CMAKE_ARGS "${google_benchmark_args}"
    UPDATE_DISCONNECTED yes
)

add_library (benchmark STATIC IMPORTED)
add_dependencies(benchmark google_benchmark)
set_target_properties(benchmark PROPERTIES IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/lib/libbenchmark.a")

macro(seqan3_benchmark target_cpp)
    # $target_cpp = "pod_tuple_benchmark.cpp"
    #   * will register the global TARGET name "pod_tuple_benchmark" and
    #   * will register the benchmark name "core_pod_tuple" if
    #     pod_tuple_benchmark.cpp is in benchmark/core/
    #
    # NOTE(marehr): ".+/benchmark/" REGEX is greedy, that means
    # /benchmark/benchmark/benchmark/hello_benchmark.cpp will result in an empty `benchmark_path`
    string(REGEX REPLACE "_benchmark.cpp$" "" target_name ${target_cpp})
    string(REGEX REPLACE ".+/benchmark/" "" benchmark_path ${CMAKE_CURRENT_LIST_DIR})
    string(REGEX REPLACE "/" "_" benchmark_path ${benchmark_path})

    add_executable(${target_name}_benchmark ${target_cpp})
    target_link_libraries(${target_name}_benchmark "pthread" "benchmark")
    add_dependencies(${target_name}_benchmark benchmark)
endmacro()

# ----------------------------------------------------------------------------
# Add the benchmarks
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
