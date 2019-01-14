# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
# -----------------------------------------------------------------------------------------------------

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
    set_target_properties(${target}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${source_path}"
    )

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
