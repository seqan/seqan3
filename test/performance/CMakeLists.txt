# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.2)
project (seqan3_test_performance CXX)

include (../seqan3-test.cmake)

macro (seqan3_benchmark target_cpp)
    # $target_cpp = "pod_tuple_benchmark.cpp"
    #   * will register the global TARGET name "pod_tuple_benchmark" and
    #   * will register the benchmark name "core_pod_tuple" if
    #     pod_tuple_benchmark.cpp is in benchmark/core/
    #
    # NOTE(marehr): ".+/benchmark/" REGEX is greedy, that means
    # /benchmark/benchmark/benchmark/hello_benchmark.cpp will result in an empty `benchmark_path`
    string (REGEX REPLACE "_benchmark.cpp$" "" target_name ${target_cpp})
    string (REGEX REPLACE ".+/test/" "" benchmark_path ${CMAKE_CURRENT_LIST_DIR})
    set (target "${target_name}_benchmark")

    add_executable (${target} ${target_cpp})
    target_link_libraries (${target} seqan3::test::performance)
    add_test (NAME "${benchmark_path}/${target}" COMMAND ${target})

    unset (target_name)
    unset (benchmark_path)
    unset (target)
endmacro ()

seqan3_require_ccache ()
seqan3_require_benchmark ()

add_subdirectories ()
