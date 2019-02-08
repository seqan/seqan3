# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.2)
project (seqan3_test_unit CXX)

include (../seqan3-test.cmake)

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
    set (target "${target_name}_test")

    add_executable (${target} ${target_cpp})
    target_link_libraries (${target} seqan3::test::unit)
    add_test (NAME "${test_path}/${target}" COMMAND ${target})

    unset (target_name)
    unset (test_path)
    unset (target)
endmacro ()

seqan3_require_ccache ()
seqan3_require_test ()

add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/test)
