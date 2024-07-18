# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# This file provides functionality common to the different test modules used by
# SeqAn3. To build tests, run cmake on one of the sub-folders in this directory
# which contain a CMakeLists.txt.

cmake_minimum_required (VERSION 3.10)

# require SeqAn3 package
find_package (SeqAn3 REQUIRED HINTS ${CMAKE_CURRENT_LIST_DIR}/../build_system)

include (CheckCXXCompilerFlag)
include (CheckCXXSourceCompiles)
include (FindPackageHandleStandardArgs)
include (FindPackageMessage)

option (SEQAN3_TEST_BUILD_OFFLINE "Skip the update step of external projects." OFF)

# Force alignment of benchmarked loops so that numbers are reliable.
# For large loops and erratic seeming bench results the value might
# have to be adapted or the option deactivated.
option (SEQAN3_BENCHMARK_ALIGN_LOOPS "Pass -falign-loops=32 to the benchmark builds." ON)

# ----------------------------------------------------------------------------
# Custom Build types
# ----------------------------------------------------------------------------

# -DCMAKE_BUILD_TYPE=FEDORA; our library did not compile for fedora quite a few times, because of that we created this
# custom build type to emulate their flag set-up.
# We omitted:
#   -specs=/usr/lib/rpm/redhat/redhat-hardened-cc1
#   -specs=/usr/lib/rpm/redhat/redhat-annobin-cc1
#   set -mtune=native
#   and -fcf-protection=check
# See https://src.fedoraproject.org/rpms/redhat-rpm-config/blob/rawhide/f/buildflags.md for an overview
set (CMAKE_CXX_FLAGS_FEDORA
     "-O2 -flto -ffat-lto-objects -fexceptions -g -grecord-gcc-switches -pipe -Wall -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -Wp,-D_GLIBCXX_ASSERTIONS -fstack-protector-strong -m64 -mtune=native -fasynchronous-unwind-tables -fstack-clash-protection -fcf-protection=check"
)

# ----------------------------------------------------------------------------
# Paths to folders.
# ----------------------------------------------------------------------------

find_path (SEQAN3_TEST_INCLUDE_DIR
           NAMES seqan3/test/tmp_directory.hpp
           HINTS "${CMAKE_CURRENT_LIST_DIR}/include/")
find_path (SEQAN3_TEST_CMAKE_MODULE_DIR
           NAMES seqan3_test_component.cmake
           HINTS "${CMAKE_CURRENT_LIST_DIR}/cmake/")
list (APPEND CMAKE_MODULE_PATH "${SEQAN3_TEST_CMAKE_MODULE_DIR}")

# ----------------------------------------------------------------------------
# Interface targets for the different test modules in seqan3.
# ----------------------------------------------------------------------------

# seqan3::test exposes a base set of required flags, includes, definitions and
# libraries which are in common for **all** seqan3 tests
if (NOT TARGET seqan3::test)
    add_library (seqan3_test INTERFACE)
    target_compile_options (seqan3_test INTERFACE "-pedantic" "-Wall" "-Wextra" "-Werror")

    # GCC12 and above: Disable warning about std::hardware_destructive_interference_size not being ABI-stable.
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12)
            target_compile_options (seqan3_test INTERFACE "-Wno-interference-size")
        endif ()
    endif ()

    target_link_libraries (seqan3_test INTERFACE "seqan3::seqan3" "pthread")
    target_include_directories (seqan3_test INTERFACE "${SEQAN3_TEST_INCLUDE_DIR}")
    add_library (seqan3::test ALIAS seqan3_test)
endif ()

# seqan3::test::performance specifies required flags, includes and libraries
# needed for performance test cases in seqan3/test/performance
if (NOT TARGET seqan3::test::performance)
    add_library (seqan3_test_performance INTERFACE)
    target_link_libraries (seqan3_test_performance INTERFACE "seqan3::test" "benchmark_main" "benchmark")
    # std::views::join is experimental in libc++
    target_compile_definitions (seqan3_test_performance INTERFACE _LIBCPP_ENABLE_EXPERIMENTAL)

    check_cxx_compiler_flag ("-falign-loops=32" SEQAN3_HAS_FALIGN_LOOPS)
    if (SEQAN3_BENCHMARK_ALIGN_LOOPS AND SEQAN3_HAS_FALIGN_LOOPS)
        target_compile_options (seqan3_test_performance INTERFACE "-falign-loops=32")
    endif ()

    add_library (seqan3::test::performance ALIAS seqan3_test_performance)
endif ()

# seqan3::test::unit specifies required flags, includes and libraries
# needed for unit test cases in seqan3/test/unit
if (NOT TARGET seqan3::test::unit)
    add_library (seqan3_test_unit INTERFACE)
    target_link_libraries (seqan3_test_unit INTERFACE "seqan3::test" "gtest_main" "gtest")
    add_library (seqan3::test::unit ALIAS seqan3_test_unit)
endif ()

# seqan3::test::header specifies required flags, includes and libraries
# needed for header test cases in seqan3/test/header
if (NOT TARGET seqan3::test::header)
    add_library (seqan3_test_header INTERFACE)
    target_link_libraries (seqan3_test_header INTERFACE "seqan3::test::unit")
    target_link_libraries (seqan3_test_header INTERFACE "seqan3::test::performance")
    target_compile_options (seqan3_test_header INTERFACE "-Wno-unused-const-variable" "-Wno-unused-variable"
                                                         "-Wno-unused-function")
    target_compile_definitions (seqan3_test_header INTERFACE -DSEQAN3_DISABLE_DEPRECATED_WARNINGS)
    target_compile_definitions (seqan3_test_header INTERFACE -DSEQAN3_HEADER_TEST)
    add_library (seqan3::test::header ALIAS seqan3_test_header)
endif ()

# ----------------------------------------------------------------------------
# Commonly shared options for external projects.
# ----------------------------------------------------------------------------

set (SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")

# ----------------------------------------------------------------------------
# Commonly used macros for the different test modules in seqan3.
# ----------------------------------------------------------------------------

include (seqan3_test_component)
include (seqan3_test_files)
include (seqan3_require_ccache)
include (seqan3_require_benchmark)
include (seqan3_require_test)
include (add_subdirectories)
