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
project (seqan3_header_test CXX)

include (../seqan3-test.cmake)

# We compile each header twice in separate compilation units. Each alone is
# sufficient to test that the header is functional, but both are needed to check
# for link errors, which can happen if the header accidentally defines a
# variable, e.g. a global or class static member. Furthermore this tests that
# header guards are working by including the same header twice.
#
# example invocation:
#     seqan3_header_test (seqan3 "<path>/include/seqan3")
#
# \sa Modified version from Bio-Formats
# https://github.com/openmicroscopy/bioformats/blob/d3bb33eeda23e81f78fd25f658bfc14a4363805f/cpp/cmake/HeaderTest.cmake#L81-L113
macro (seqan3_header_test component header_base_path)
    set (target "${component}_header_test")
    separate_arguments (cxx_flags_list UNIX_COMMAND "${SEQAN3_CXX_FLAGS}")

    file (RELATIVE_PATH header_relative_base_path "${CMAKE_SOURCE_DIR}" "${header_base_path}")

    # finding all *.hpp files relative from the current directory (e.g. /test/)
    # The resulting list is normalized to `header_base_path` that means concatenating
    # "${header_base_path}/header_files[i]" will result in an absolute path to the file
    #
    # Example output:
    #   seqan3/alphabet/adaptation/all.hpp
    #   seqan3/alphabet/adaptation/char.hpp
    #   seqan3/alphabet/adaptation/concept.hpp
    #   seqan3/alphabet/adaptation/uint.hpp
    #   seqan3/alphabet/all.hpp
    #   seqan3/alphabet/dna5_detail.hpp <- will be filtered out
    #   ....
    file (GLOB_RECURSE header_files RELATIVE "${header_base_path}"
         "${header_relative_base_path}/*.hpp"
         "${header_relative_base_path}/*.h")

    # filter out headers which are implementation detail
    list (FILTER header_files EXCLUDE REGEX "detail|/bits/")

    file (MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/header_test_src")

    file (WRITE "${PROJECT_BINARY_DIR}/header_test_src/empty.cpp" "")

    add_executable (${target} ${PROJECT_BINARY_DIR}/header_test_src/empty.cpp)
    target_compile_options (${target} PRIVATE ${cxx_flags_list} ${SEQAN3_TEST_CXX_FLAGS})
    target_compile_definitions (${target} PRIVATE ${SEQAN3_DEFINITIONS} ${SEQAN3_TEST_DEFINITIONS})
    target_include_directories (${target} PRIVATE ${SEQAN3_INCLUDE_DIRS} ${SEQAN3_TEST_INCLUDE_DIRS})
    target_link_libraries (${target} ${SEQAN3_LIBRARIES} ${SEQAN3_TEST_LIBRARIES})
    add_test (NAME "header/${target}" COMMAND ${target})
    add_dependencies (${target} gtest)

    foreach (header ${header_files})
        foreach (repeat 1 2)
            # header = seqan3/alphabet/adaptation/all.hpp
            # header_source_file_name = seqan3_alphabet_adaptation_all_1.cpp
            # header_source_test_name = seqan3_alphabet_adaptation_all_1
            string (REPLACE "/" "_" header_source_file_name ${header})
            string (REGEX REPLACE "\\.(hpp|h)$" "_${repeat}.cpp" header_source_file_name ${header_source_file_name})
            string (REGEX REPLACE "\\.cpp$" "" header_source_test_name ${header_source_file_name})
            set (target_library "${component}_header_${header_source_test_name}")
            set (source_file "${PROJECT_BINARY_DIR}/header_test_src/${header_source_file_name}")

            file (WRITE "${source_file}" "
#include <gtest/gtest.h>
#include <${header}>
#include <${header}>
TEST(Header, ${header_source_test_name}) {}")

            add_library (${target_library} OBJECT "${source_file}")
            target_compile_options (${target_library} PRIVATE ${cxx_flags_list} ${SEQAN3_TEST_CXX_FLAGS})
            target_compile_definitions (${target_library} PRIVATE ${SEQAN3_DEFINITIONS} ${SEQAN3_TEST_DEFINITIONS})
            target_include_directories (${target_library} PRIVATE ${SEQAN3_INCLUDE_DIRS} ${SEQAN3_TEST_INCLUDE_DIRS})
            add_dependencies (${target_library} gtest)

            target_sources (${target} PRIVATE $<TARGET_OBJECTS:${target_library}>)
        endforeach ()
    endforeach ()

    unset (target)
    unset (cxx_flags_list)
    unset (header_relative_base_path)
    unset (header_files)
    unset (header_source_file_name)
    unset (header_source_test_name)
    unset (target_library)
    unset (source_file)
endmacro ()

seqan3_require_ccache ()
seqan3_require_test ()

seqan3_header_test (seqan3 "${SEQAN3_CLONE_DIR}/include")
seqan3_header_test (seqan3_test "${SEQAN3_CLONE_DIR}/test/include")
seqan3_header_test (range-v3 "${SEQAN3_CLONE_DIR}/submodules/range-v3/include")

# not self-contained headers; error: extra ‘;’ [-Werror=pedantic]
# seqan3_header_test (lemon "${SEQAN3_CLONE_DIR}/submodules/lemon/include")

# not self-contained headers
# seqan3_header_test (cereal "${SEQAN3_CLONE_DIR}/submodules/cereal/include")

# not self-contained headers; multiple defined symbols
# seqan3_header_test (sdsl-lite "${SEQAN3_CLONE_DIR}/submodules/sdsl-lite/include")
