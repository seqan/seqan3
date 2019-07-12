# Minimum cmake version
cmake_minimum_required(VERSION 3.7)

# Name of project
project(seqan3_documentation NONE)

### Find seqan3
if (NOT SEQAN3_INCLUDE_DIR)
    find_path (SEQAN3_INCLUDE_DIR NAMES seqan3/version.hpp HINTS "${CMAKE_SOURCE_DIR}/../../include")
endif ()

if (NOT EXISTS "${SEQAN3_INCLUDE_DIR}/seqan3/version.hpp")
    message (FATAL_ERROR "Could not find SeqAn3. Not building documentation.")
endif ()
message (STATUS "Found: ${SEQAN3_INCLUDE_DIR}")

set (SEQAN3_DOXYGEN_INPUT_DIR "${CMAKE_SOURCE_DIR}")

include (seqan3-doxygen.cmake)

### Find python interpreter
find_package(PythonInterp REQUIRED)

### Enable testing
enable_testing()

set (SEQAN3_DOXYGEN_TEST_SCRIPT "
    output=\$(${DOXYGEN_EXECUTABLE} 2>&1 > /dev/null | sort);
    echo \"\$output\";
    test -z \"\$output\"")

if (SEQAN3_USER_DOC)
    message (STATUS "Add user doc test.")
    add_test(NAME doc_usr_test
    COMMAND sh -c "${SEQAN3_DOXYGEN_TEST_SCRIPT}"
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/doc_usr)
endif ()

if (SEQAN3_DEV_DOC)
    message (STATUS "Add dev doc test.")
    add_test(NAME doc_dev_test
    COMMAND sh -c "${SEQAN3_DOXYGEN_TEST_SCRIPT}"
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/doc_dev)
endif ()
