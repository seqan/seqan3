# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

# Minimum cmake version
cmake_minimum_required(VERSION 3.7)

### Find seqan3
include (../../build_system/seqan3-config-version.cmake)
set (SEQAN3_VERSION "${PACKAGE_VERSION}")

project (seqan3 LANGUAGES NONE VERSION "${SEQAN3_VERSION}")

if (NOT EXISTS "${SEQAN3_INCLUDE_DIR}/seqan3/version.hpp")
    message (FATAL_ERROR "Could not find SeqAn3. Not building documentation.")
endif ()
message (STATUS "Found: ${SEQAN3_INCLUDE_DIR} (${SEQAN3_VERSION})")

set (SEQAN3_DOXYGEN_INPUT_DIR "${CMAKE_SOURCE_DIR}")

include (seqan3-doxygen.cmake)

enable_testing ()

if (SEQAN3_USER_DOC)
    add_subdirectory(doc_usr)
endif ()

if (SEQAN3_DEV_DOC)
    add_subdirectory(doc_dev)
endif ()

include (seqan3-doxygen-package.cmake)
