# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

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

if (SEQAN3_USER_DOC)
    add_subdirectory (doc_usr)
endif ()

if (SEQAN3_DEV_DOC)
    add_subdirectory (doc_dev)
endif ()

include (seqan3-doxygen-package.cmake)
