# Minimum cmake version
cmake_minimum_required(VERSION 3.7)

# Name of project
project(seqan3_documentation NONE)

### Find seqan3
include (../../build_system/seqan3-config-version.cmake)
set (SEQAN3_VERSION "${PACKAGE_VERSION}")

if (NOT EXISTS "${SEQAN3_INCLUDE_DIR}/seqan3/version.hpp")
    message (FATAL_ERROR "Could not find SeqAn3. Not building documentation.")
endif ()
message (STATUS "Found: ${SEQAN3_INCLUDE_DIR} (${SEQAN3_VERSION})")

set (SEQAN3_DOXYGEN_INPUT_DIR "${CMAKE_SOURCE_DIR}")

include (seqan3-doxygen.cmake)

### Enable testing
enable_testing()

# doxygen does not show any warnings (doxygen prints warnings / errors to cerr)
set (SEQAN3_DOXYGEN_FAIL_ON_WARNINGS "
    ${DOXYGEN_EXECUTABLE} > doxygen.cout 2> doxygen.cerr;
    cat \"doxygen.cerr\";
    test ! -s \"doxygen.cerr\"")

# We search the HTML output to ensure that no `requires` clauses are at wrong places.
set (SEQAN3_DOXYGEN_FAIL_ON_UNCOND_REQUIRES
     "! find . -not -name \"*_source.html\" -name \"*.html\" -print0 | xargs -0 grep \"requires\" | grep \"memname\"")

if (SEQAN3_USER_DOC)
    message (STATUS "Add user doc test.")
    add_test(NAME doc_usr_no_warnings_test
             COMMAND bash -c "${SEQAN3_DOXYGEN_FAIL_ON_WARNINGS}"
             WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/doc_usr)

    add_test(NAME doc_usr_no_uncond_requires_test
             COMMAND bash -c "${SEQAN3_DOXYGEN_FAIL_ON_UNCOND_REQUIRES}"
             WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/doc_usr)
endif ()

if (SEQAN3_DEV_DOC)
    message (STATUS "Add dev doc test.")
    add_test(NAME doc_dev_no_warnings_test
             COMMAND bash -c "${SEQAN3_DOXYGEN_FAIL_ON_WARNINGS}"
             WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/doc_dev)

    add_test(NAME doc_dev_no_uncond_requires_test
             COMMAND bash -c "${SEQAN3_DOXYGEN_FAIL_ON_UNCOND_REQUIRES}"
             WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/doc_dev)
endif ()
