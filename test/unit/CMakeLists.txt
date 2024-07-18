# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.10...3.22)
project (seqan3_test_unit CXX)

include (../seqan3-test.cmake)
include (GoogleTest OPTIONAL)

include (diagnostics/list_missing_unit_tests)
include (diagnostics/list_unused_unit_tests)
include (include_dependencies/add_include_dependencies)

option (SEQAN3_VERBOSE_TESTS "Run each test case individually" OFF)
option (SEQAN3_USE_INCLUDE_DEPENDENCIES
        "Build tests in an hierarchical order (by an include graph, i.e. tests with less dependencies are build first)"
        OFF)

macro (seqan3_test unit_test_cpp)
    cmake_parse_arguments (SEQAN3_TEST "" "" "CYCLIC_DEPENDING_INCLUDES" ${ARGN})

    file (RELATIVE_PATH unit_test "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${unit_test_cpp}")
    seqan3_test_component (target "${unit_test}" TARGET_NAME)
    seqan3_test_component (test_name "${unit_test}" TEST_NAME)

    add_executable (${target} ${unit_test_cpp})
    target_link_libraries (${target} seqan3::test::unit)
    add_include_dependencies (${target} "${SEQAN3_TEST_CYCLIC_DEPENDING_INCLUDES}")
    collect_used_test (${target})
    if (SEQAN3_VERBOSE_TESTS)
        gtest_discover_tests (
            ${target}
            TEST_PREFIX "${test_name}::"
            PROPERTIES
            TIMEOUT "30")
    else ()
        add_test (NAME "${test_name}" COMMAND ${target})
    endif ()

    unset (unit_test)
    unset (target)
    unset (test_name)
endmacro ()

seqan3_require_ccache ()
seqan3_require_test ()

add_subdirectories ()

list_missing_unit_tests ()
list_unused_unit_tests ()
