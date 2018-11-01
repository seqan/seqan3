# Minimum cmake version
cmake_minimum_required(VERSION 3.0.0)

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

if (SEQAN3_USER_DOC)
    message (STATUS "Add user doc test.")
    add_test(NAME test_user_doc
             COMMAND ${PYTHON_EXECUTABLE}
             ${CMAKE_SOURCE_DIR}/doc_test.py "-e" ${DOXYGEN_EXECUTABLE} "-i" ${seqan3_doxyfile_user})
endif ()

if (SEQAN3_DEV_DOC)
    message (STATUS "Add dev doc test.")
    add_test(NAME test_devel_doc
             COMMAND ${PYTHON_EXECUTABLE}
             ${CMAKE_SOURCE_DIR}/doc_test.py "-e" ${DOXYGEN_EXECUTABLE} "-i" ${seqan3_doxyfile_devel})
endif ()
