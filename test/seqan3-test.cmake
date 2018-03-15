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

# require SeqAn3 package
find_package (SeqAn3 REQUIRED
              HINTS ${CMAKE_CURRENT_LIST_DIR}/../build_system)

include(CheckCXXSourceCompiles)
include(FindPackageHandleStandardArgs)
include(FindPackageMessage)

set(SEQAN3_BENCHMARK_SRC_DIR "${PROJECT_BINARY_DIR}/vendor/benchmark")
set(SEQAN3_TEST_SRC_DIR "${PROJECT_BINARY_DIR}/vendor/googletest")

# required flags, includes, definitions and libraries for seqan3
set(SEQAN3_STRICT_CXX_FLAGS "-pedantic" "-Werror" "-Wall" "-Wextra")

# required flags, includes and libraries for seqan3/test/performance
set(SEQAN3_BENCHMARK_CXX_FLAGS ${SEQAN3_STRICT_CXX_FLAGS})
set(SEQAN3_BENCHMARK_INCLUDE_DIRS "")
set(SEQAN3_BENCHMARK_LIBRARIES "")
list(APPEND SEQAN3_BENCHMARK_INCLUDE_DIRS "${SEQAN3_BENCHMARK_SRC_DIR}/include/")
list(APPEND SEQAN3_BENCHMARK_INCLUDE_DIRS "${SEQAN3_SRC_DIR}/test/include/")

# required flags, includes and libraries for seqan3/test/unit
set(SEQAN3_TEST_CXX_FLAGS ${SEQAN3_STRICT_CXX_FLAGS})
set(SEQAN3_TEST_INCLUDE_DIRS "")
set(SEQAN3_TEST_LIBRARIES "")
list(APPEND SEQAN3_TEST_INCLUDE_DIRS "${SEQAN3_TEST_SRC_DIR}/googletest/include/")
list(APPEND SEQAN3_TEST_INCLUDE_DIRS "${SEQAN3_SRC_DIR}/test/include/")

# commonly shared options:
set(SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "")
list(APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list(APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list(APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")

macro(seqan3_require_ccache)
    find_program(CCACHE_PROGRAM ccache)
    find_package_message(CCACHE_PROGRAM_PRE "Finding program ccache" "[${CCACHE_PROGRAM}]")

    if(NOT CCACHE_PROGRAM)
        find_package_message(CCACHE_PROGRAM "Finding program ccache - Failed" "[${CCACHE_PROGRAM}]")
    else()
        find_package_message(CCACHE_PROGRAM "Finding program ccache - Success" "[${CCACHE_PROGRAM}]")
        if(CMAKE_VERSION VERSION_LESS 3.4)
            set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
            set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
        else()
            # New option since cmake >= 3.4:
            # https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_LAUNCHER.html
            set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")

            # use ccache in external cmake projects
            list(APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER}")
        endif()
    endif()
    unset(CCACHE_PROGRAM)
endmacro()

macro(seqan3_require_benchmark)
    enable_testing()

    set(google_benchmark_args ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS})
    list(APPEND google_benchmark_args "-DBENCHMARK_ENABLE_TESTING=false")
    # list(APPEND google_benchmark_args "-DBENCHMARK_ENABLE_LTO=true")

    include(ExternalProject)
    ExternalProject_Add(
        google_benchmark
        PREFIX google_benchmark
        GIT_REPOSITORY "https://github.com/google/benchmark.git"
        SOURCE_DIR "${SEQAN3_BENCHMARK_SRC_DIR}"
        CMAKE_ARGS "${google_benchmark_args}"
        UPDATE_DISCONNECTED yes
    )
    unset(google_benchmark_args)

    set(benchmark_library "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}benchmark${CMAKE_STATIC_LIBRARY_SUFFIX}")

    add_library (benchmark STATIC IMPORTED)
    add_dependencies(benchmark google_benchmark)
    set_target_properties(benchmark PROPERTIES IMPORTED_LOCATION "${benchmark_library}")

    list(APPEND SEQAN3_BENCHMARK_LIBRARIES "benchmark")
    list(APPEND SEQAN3_BENCHMARK_LIBRARIES "pthread")
endmacro()

macro(seqan3_require_test)
    enable_testing()

    set(google_test_args ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS})
    list(APPEND google_test_args "-DBUILD_GTEST=1")
    list(APPEND google_test_args "-DBUILD_GMOCK=0")

    # force that libraries are installed to `lib/`, because GNUInstallDirs might install it into `lib64/`
    list(APPEND google_test_args "-DCMAKE_INSTALL_LIBDIR=${PROJECT_BINARY_DIR}/lib/")

    include(ExternalProject)
    ExternalProject_Add(
        googletest
        PREFIX googletest
        GIT_REPOSITORY "https://github.com/google/googletest.git"
        GIT_TAG "15392f1a38fa0b8c3f13a9732e94b209069efa1c"
        SOURCE_DIR "${SEQAN3_TEST_SRC_DIR}"
        CMAKE_ARGS "${google_test_args}"
        UPDATE_DISCONNECTED yes
    )
    unset(google_test_args)

    # google sets CMAKE_DEBUG_POSTFIX = "d"
    set(gtest_main_library "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX}")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(gtest_main_library "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_maind${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif()

    add_library (gtest_main STATIC IMPORTED)
    add_dependencies(gtest_main googletest)
    set_target_properties(gtest_main PROPERTIES IMPORTED_LOCATION "${gtest_main_library}")

    set(gtest_library "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(gtest_library "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtestd${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif()

    add_library (gtest STATIC IMPORTED)
    add_dependencies(gtest gtest_main)
    set_target_properties(gtest PROPERTIES IMPORTED_LOCATION "${gtest_library}")

    list(APPEND SEQAN3_TEST_LIBRARIES "gtest_main")
    list(APPEND SEQAN3_TEST_LIBRARIES "gtest")
    list(APPEND SEQAN3_TEST_LIBRARIES "pthread")
endmacro()

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
    unset(ENTRIES)
endmacro ()
