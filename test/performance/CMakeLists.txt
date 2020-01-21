# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.7)
project (seqan3_test_performance CXX)

include (../seqan3-test.cmake)

macro (seqan3_benchmark benchmark_cpp)
    file (RELATIVE_PATH benchmark "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${benchmark_cpp}")
    seqan3_test_component (target "${benchmark}" TARGET_NAME)
    seqan3_test_component (test_name "${benchmark}" TEST_NAME)

    # Force alignment of benchmarked loops so that numbers are reliable.
    # For large loops and erratic seeming bench results the value might
    # have to be adapted or the option deactivated.
    option (BENCHMARK_ALIGN_LOOPS "Pass -falign-loops=32 to the benchmark builds." ON)

    add_executable (${target} ${benchmark_cpp})
    target_link_libraries (${target} seqan3::test::performance)
    if (BENCHMARK_ALIGN_LOOPS)
        target_compile_options (${target} PUBLIC "-falign-loops=32")
    endif ()
    add_test (NAME "${test_name}" COMMAND ${target})

    unset (benchmark)
    unset (target)
    unset (test_name)
endmacro ()

seqan3_require_ccache ()
seqan3_require_benchmark ()

add_subdirectories ()
