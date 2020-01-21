# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.7)
project (seqan3_test_coverage CXX)

include (../seqan3-test.cmake)

if (CMAKE_BUILD_TYPE AND NOT ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug"))
    message(WARNING "Coverage test must be build in debug mode [build type = ${CMAKE_BUILD_TYPE}]")
endif ()

find_program (LCOV_COMMAND  NAMES lcov lcov.bat lcov.exe lcov.perl)
find_program (GENHTML_COMMAND NAMES genhtml genhtml.perl genhtml.bat)

if (NOT LCOV_COMMAND)
    message(FATAL_ERROR "lcov not found! Aborting...")
endif ()

# Files which should not be included in the coverage report
set (TEST_COVERAGE_EXCLUDE_FILES
    "'/usr/*'"
    "'${SEQAN3_CLONE_DIR}/include/seqan3/contrib/*'"
    "'${SEQAN3_CLONE_DIR}/include/seqan3/std/*'"
    "'${SEQAN3_CLONE_DIR}/submodules/*'"
    "'${SEQAN3_CLONE_DIR}/test/unit/*'"
    "'${PROJECT_BINARY_DIR}/vendor/*'"
)
# Holds all target's defined by seqan3_test
set_property (GLOBAL PROPERTY GLOBAL_TEST_COVERAGE_ALL_TESTS "")

add_custom_command (
    OUTPUT ${PROJECT_BINARY_DIR}/seqan3_coverage
    # Cleanup lcov (resetting code coverage counters to zero)
    COMMAND ${LCOV_COMMAND} --directory ${PROJECT_BINARY_DIR} --zerocounters
    # Create baseline to make sure untouched files show up in the report
    COMMAND ${LCOV_COMMAND} --directory ${PROJECT_BINARY_DIR} --capture --initial --output-file ${PROJECT_BINARY_DIR}/seqan3_coverage.baseline

    # Run tests
    COMMAND ${CMAKE_CTEST_COMMAND}

    # Capturing lcov counters and generating report
    COMMAND ${LCOV_COMMAND} --directory ${PROJECT_BINARY_DIR} --capture --output-file ${PROJECT_BINARY_DIR}/seqan3_coverage.captured
    # merge baseline counters and captured counters
    COMMAND ${LCOV_COMMAND} -a ${PROJECT_BINARY_DIR}/seqan3_coverage.baseline -a ${PROJECT_BINARY_DIR}/seqan3_coverage.captured --output-file ${PROJECT_BINARY_DIR}/seqan3_coverage.total
    COMMAND ${LCOV_COMMAND} --remove ${PROJECT_BINARY_DIR}/seqan3_coverage.total ${TEST_COVERAGE_EXCLUDE_FILES} --output-file ${PROJECT_BINARY_DIR}/seqan3_coverage

    BYPRODUCTS seqan3_coverage.baseline seqan3_coverage.captured seqan3_coverage.total

    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Processing code coverage counters."
)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "seqan3_coverage.baseline;seqan3_coverage.captured;seqan3_coverage.total")

add_custom_target (
    coverage ALL
    DEPENDS
    ${PROJECT_BINARY_DIR}/seqan3_coverage
)

add_custom_target (
    coverage_html
    COMMAND ${GENHTML_COMMAND} --highlight --legend --output-directory html ${PROJECT_BINARY_DIR}/seqan3_coverage
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    DEPENDS coverage
    COMMENT "Generate coverage report."
)

add_custom_command(
    TARGET coverage_html POST_BUILD
    COMMAND ;
    COMMENT "Open ${PROJECT_BINARY_DIR}/html/index.html in your browser to view the coverage report."
)

macro (seqan3_test unit_test_cpp)
    file (RELATIVE_PATH unit_test "${CMAKE_SOURCE_DIR}/../unit" "${CMAKE_CURRENT_LIST_DIR}/${unit_test_cpp}")
    seqan3_test_component (target "${unit_test}" TARGET_NAME)
    seqan3_test_component (test_name "${unit_test}" TEST_NAME)

    add_executable (${target} ${unit_test_cpp})
    target_link_libraries (${target} seqan3::test::coverage)
    add_test (NAME "${test_name}" COMMAND ${target})

    # any change of a target will invalidate the coverage result;
    # NOTE that this is a GLOBAL variable, because a normal
    # `set(GLOBAL_TEST_COVERAGE_ALL_TESTS)` would not propagate the result when
    # CMakeLists.txt goes out of scope due to a `add_subdirectory`
    set_property(GLOBAL APPEND PROPERTY GLOBAL_TEST_COVERAGE_ALL_TESTS ${target})

    unset (unit_test)
    unset (target)
    unset (test_name)
endmacro ()

seqan3_require_ccache ()
seqan3_require_test ()

# add all unit tests
add_subdirectories_of ("${CMAKE_CURRENT_SOURCE_DIR}/../unit")

# add collected test cases as dependency
get_property(TEST_COVERAGE_ALL_TESTS GLOBAL PROPERTY GLOBAL_TEST_COVERAGE_ALL_TESTS)
add_custom_command (
    OUTPUT ${PROJECT_BINARY_DIR}/seqan3_coverage
    DEPENDS ${TEST_COVERAGE_ALL_TESTS}
    APPEND
)
