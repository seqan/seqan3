# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.20...3.31)
project (seqan3_test_coverage CXX)

# Add a custom build type: Coverage
#
# `--coverage` is equivalent to `-fprofile-arcs -ftest-coverage`
# https://gcc.gnu.org/onlinedocs/gcc/Instrumentation-Options.html#index-gcov
#
# `-fprofile-update=atomic` prevents profile corruption in multi-threaded tests
# https://gcc.gnu.org/onlinedocs/gcc/Instrumentation-Options.html#index-fprofile-update
#
# `-fprofile-abs-path` converts relative source file names to absolute paths in the coverage files
# https://gcc.gnu.org/onlinedocs/gcc/Instrumentation-Options.html#index-fprofile-abs-path

set (CMAKE_CXX_FLAGS_COVERAGE
     "${CMAKE_CXX_FLAGS_DEBUG} --coverage -g -O0 -fprofile-abs-path -fprofile-update=atomic"
     CACHE STRING "Flags used by the C++ compiler during coverage builds." FORCE)
set (CMAKE_C_FLAGS_COVERAGE
     "${CMAKE_C_FLAGS_DEBUG} --coverage -g -O0 -fprofile-abs-path -fprofile-update=atomic"
     CACHE STRING "Flags used by the C compiler during coverage builds." FORCE)
set (CMAKE_EXE_LINKER_FLAGS_COVERAGE
     "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -Wl,-lgcov"
     CACHE STRING "Flags used for linking binaries during coverage builds." FORCE)
set (CMAKE_SHARED_LINKER_FLAGS_COVERAGE
     "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -Wl,-lgcov"
     CACHE STRING "Flags used by the shared libraries linker during coverage builds." FORCE)

mark_as_advanced (CMAKE_CXX_FLAGS_COVERAGE CMAKE_C_FLAGS_COVERAGE CMAKE_EXE_LINKER_FLAGS_COVERAGE
                  CMAKE_SHARED_LINKER_FLAGS_COVERAGE)

enable_testing ()
add_subdirectory ("../unit" "${CMAKE_CURRENT_BINARY_DIR}/unit")
