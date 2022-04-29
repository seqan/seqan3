# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

# Minimum cmake version
cmake_minimum_required (VERSION 3.7)

### Find seqan3
include (../../build_system/seqan3-config-version.cmake)
set (SEQAN3_VERSION "${PACKAGE_VERSION}")

project (seqan3
         LANGUAGES NONE
         VERSION "${SEQAN3_PROJECT_VERSION}")

if (NOT EXISTS "${SEQAN3_INCLUDE_DIR}/seqan3/version.hpp")
    message (FATAL_ERROR "Could not find SeqAn3. Not building documentation.")
endif ()
message (STATUS "Found: ${SEQAN3_INCLUDE_DIR} (${SEQAN3_VERSION})")

set (SEQAN3_DOXYGEN_INPUT_DIR "${CMAKE_SOURCE_DIR}")

include (seqan3-doxygen.cmake)

include (seqan3-doxygen-layout.cmake)

enable_testing ()

### Copy over snippets from sharg
### -----------------------------
### This is needed to replace sharg with seqan3 in all snippets for correct documentation.
set (SHARG_SNIPPET_OUTPUT_DIR "${CMAKE_BINARY_DIR}/sharg_snippets")
file(MAKE_DIRECTORY ${SHARG_SNIPPET_OUTPUT_DIR})
message(STATUS "Created sharg snippet directory: ${SHARG_SNIPPET_OUTPUT_DIR}")

include (${CMAKE_SOURCE_DIR}/../cmake/seqan3_test_files.cmake)
seqan3_test_files (sharg_snippet_files "${SEQAN3_INCLUDE_DIR}/../submodules/sharg-parser/test/snippet/" "*.cpp")

foreach (sharg_snippet_filename ${sharg_snippet_files})
    file(READ "${SEQAN3_INCLUDE_DIR}/../submodules/sharg-parser/test/snippet/${sharg_snippet_filename}" SNIPPET_INPUT_CONTENT)
    string(REGEX REPLACE "sharg" "seqan3" SNIPPET_OUTPUT_CONTENT ${SNIPPET_INPUT_CONTENT})
    file(WRITE "${SHARG_SNIPPET_OUTPUT_DIR}/${sharg_snippet_filename}" ${SNIPPET_OUTPUT_CONTENT})
endforeach ()

### Add subdirectories
if (SEQAN3_USER_DOC)
    add_subdirectory (doc_usr)
endif ()

if (SEQAN3_DEV_DOC)
    add_subdirectory (doc_dev)
endif ()

include (seqan3-doxygen-package.cmake)
