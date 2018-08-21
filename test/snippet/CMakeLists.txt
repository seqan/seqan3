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
project (seqan3_test_snippet CXX)

include (../seqan3-test.cmake)

macro (seqan3_snippet target_cpp snippet_full_path)
    # $target_cpp = "pod_tuple_snippet.cpp"
    #   * will register the global TARGET name "pod_tuple_snippet" and
    #   * will register the test case name "core_pod_tuple" if
    #     pod_tuple_snippet.cpp is in test/core/
    #
    # NOTE(marehr): ".+/test/" REGEX is greedy, that means
    # /test/test/test/hello_snippet.cpp will result in an empty `test_path`
    string (REGEX REPLACE ".cpp$" "" target_name ${target_cpp})
    string (REGEX REPLACE ".+/test/" "" test_path ${snippet_full_path})
    string (REGEX REPLACE ".+/test/snippet/" "" source_path ${snippet_full_path})
    set (target "${target_name}_snippet")

    add_executable (${target} ${source_path}/${target_cpp})
    target_link_libraries (${target} seqan3::test::unit)
    # target_compile_options(${target} PUBLIC --disable-werror)

    add_test (NAME "${test_path}/${target}" COMMAND ${target})

    unset (target_name)
    unset (test_path)
    unset (target)
endmacro ()

seqan3_require_ccache ()
seqan3_require_test ()

# add all snippets
file(GLOB_RECURSE snippets "*.cpp")
foreach(snippet ${snippets})
    get_filename_component(target_cpp "${snippet}" NAME)
    get_filename_component(snippet_full_path "${snippet}" DIRECTORY)

    seqan3_snippet("${target_cpp}" "${snippet_full_path}")
endforeach()
