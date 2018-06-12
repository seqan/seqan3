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
project (seqan3_test_memory CXX)

# Load CTest module for accessing internal cmake memtest capabilities.
# Note, enable_testing() is not sufficient for this.
include (CTest)

# Load basic seqan3-test functionality
include (../seqan3-test.cmake)

# seqan3::test::header specifies required flags, includes and libraries
# needed for header test cases in seqan3/test/header
add_library (seqan3::test::memory INTERFACE IMPORTED)
set_property (TARGET seqan3::test::memory APPEND PROPERTY INTERFACE_LINK_LIBRARIES seqan3::test::unit)

# Only add memory test if valgrind is available.
find_program (valgrind_exec "valgrind")

if (valgrind_exec-NOTFOUND)
    message (FATAL_ERROR "Could not find valgrind.")
endif ()

# Redefine seqan3_test macro to add memtest.
macro (seqan3_test target_cpp)

    # $target_cpp = "pod_tuple_test.cpp"
    #   * will register the global TARGET name "pod_tuple_test" and
    #   * will register the test case name "core_pod_tuple" if
    #     pod_tuple_test.cpp is in test/core/
    #
    # NOTE(marehr): ".+/test/" REGEX is greedy, that means
    # /test/test/test/hello_test.cpp will result in an empty `test_path`
    string (REGEX REPLACE "_test.cpp$" "" target_name ${target_cpp})
    string (REGEX REPLACE ".+/test/" "" test_path ${CMAKE_CURRENT_LIST_DIR})

    # This test will generate the same binary as the standard unit test, but
    # we use it only to execute valgrind on the binary.
    set (target "${target_name}_test")

    add_executable (${target} ${target_cpp})
    target_link_libraries (${target} seqan3::test::memory)

    # We also add a custom target, that replicates the original ctest call to
    # invoke the memory debugger and print the output stored in the log.
    add_test (NAME "${test_path}/${target}" COMMAND ${target})
    set (mem_test_target "${target_name}_mem_test")
    add_custom_target(${mem_test_target}
                      COMMAND ${CMAKE_CTEST_COMMAND} --force-new-ctest-process --overwrite "MemoryCheckCommand=${valgrind_exec}" --test-action memcheck -R "${target}"
                      COMMAND cat "${CMAKE_CURRENT_BINARY_DIR}/Testing/Temporary/MemoryChecker.*.log"
                      DEPENDS "$<TARGET_FILE:${target}>"
                      COMMENT "Executing memory check ...")

    unset (target_name)
    unset (test_path)
    unset (target)
    unset (mem_test_target)
endmacro ()

seqan3_require_ccache ()
seqan3_require_test ()

add_subdirectories_of ("${CMAKE_CURRENT_SOURCE_DIR}/../unit")
