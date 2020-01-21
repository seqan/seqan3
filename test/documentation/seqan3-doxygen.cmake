# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.7)

### Find doxygen and dependency to DOT tool
message (STATUS "Searching for doxygen.")
find_package (Doxygen REQUIRED)

if (NOT ${DOXYGEN_FOUND})
    message (FATAL_ERROR "Could not find doxygen. Not building documentation.")
endif ()

if (NOT ${DOXYGEN_DOT_FOUND})
    message (STATUS "Could not find dot tool. Disabling dot support.")
    set (SEQAN3_DOXYGEN_HAVE_DOT "NO")
else ()
    message (STATUS "Found dot tool. Enabling dot support.")
    set (SEQAN3_DOXYGEN_HAVE_DOT "YES")
endif ()

### Configure doc/developer targets.
set(SEQAN3_DOXYFILE_IN ${SEQAN3_DOXYGEN_INPUT_DIR}/seqan3_doxygen_cfg.in)

option(SEQAN3_USER_DOC "Create build target and test for user documentation." ON)
option(SEQAN3_DEV_DOC "Create build target and test for developer documentation." ON)

### Download and extract cppreference-doxygen-web.tag.xml for std:: documentation links
set(SEQAN3_DOXYGEN_STD_TAGFILE "${PROJECT_BINARY_DIR}/cppreference-doxygen-web.tag.xml")
include(ExternalProject)
ExternalProject_Add (
    download-cppreference-doxygen-web-tag
    URL "https://github.com/PeterFeicht/cppreference-doc/releases/download/v20190928/html-book-20190928.tar.xz"
    URL_HASH SHA256=6c250b02f6cb0fb54ea372d81c8f50703c01cc9f6887281c88a84dd91f463abc
    TLS_VERIFY ON
    DOWNLOAD_DIR "${PROJECT_BINARY_DIR}"
    DOWNLOAD_NAME "html-book.tar.xz"
    DOWNLOAD_NO_EXTRACT YES
    BINARY_DIR "${PROJECT_BINARY_DIR}"
    CONFIGURE_COMMAND /bin/sh -c "xzcat html-book.tar.xz | tar -xf - cppreference-doxygen-web.tag.xml"
    BUILD_COMMAND rm "html-book.tar.xz"
    INSTALL_COMMAND ""
)

if (SEQAN3_USER_DOC)
    message (STATUS "Configuring user doc.")

    set (SEQAN3_DOXYGEN_OUTPUT_DIR "${PROJECT_BINARY_DIR}/doc_usr")
    set (SEQAN3_DOXYGEN_SOURCE_DIR "${SEQAN3_INCLUDE_DIR}/..")
    set (SEQAN3_DOXYGEN_EXCLUDE_SYMBOLS "detail seqan3::simd") #/""
    set (SEQAN3_DOXYGEN_PREDEFINED_NDEBUG "-NDEBUG") #/""
    set (SEQAN3_DOXYGEN_ENABLED_SECTIONS "") #/"DEV"
    set (SEQAN3_DOXYGEN_EXTRACT_PRIVATE "NO") #/"YES":

    configure_file (${SEQAN3_DOXYFILE_IN} ${SEQAN3_DOXYGEN_OUTPUT_DIR}/Doxyfile)

    add_custom_target(doc_usr ALL
                      COMMAND ${DOXYGEN_EXECUTABLE}
                      WORKING_DIRECTORY ${SEQAN3_DOXYGEN_OUTPUT_DIR}
                      DEPENDS download-cppreference-doxygen-web-tag
                      COMMENT "Generating user API documentation with Doxygen"
                      VERBATIM)
endif ()

if (SEQAN3_DEV_DOC)
    message(STATUS "Configuring devel doc.")

    set(SEQAN3_DOXYGEN_OUTPUT_DIR "${PROJECT_BINARY_DIR}/doc_dev")
    set(SEQAN3_DOXYGEN_SOURCE_DIR "${SEQAN3_INCLUDE_DIR}/..")
    set(SEQAN3_DOXYGEN_EXCLUDE_SYMBOLS "")
    set(SEQAN3_DOXYGEN_PREDEFINED_NDEBUG "")
    set(SEQAN3_DOXYGEN_ENABLED_SECTIONS "DEV")
    set(SEQAN3_DOXYGEN_EXTRACT_PRIVATE "YES")

    configure_file(${SEQAN3_DOXYFILE_IN} ${SEQAN3_DOXYGEN_OUTPUT_DIR}/Doxyfile)

    add_custom_target(doc_dev ALL
                      COMMAND ${DOXYGEN_EXECUTABLE}
                      WORKING_DIRECTORY ${SEQAN3_DOXYGEN_OUTPUT_DIR}
                      DEPENDS download-cppreference-doxygen-web-tag
                      COMMENT "Generating developer API documentation with Doxygen"
                      VERBATIM)
                      message (STATUS "Add devel doc test.")
endif ()
