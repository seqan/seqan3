# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.14)

# install and package seqan3 library
ExternalProject_Add (
    seqan3_test_prerequisite
    PREFIX seqan3_test_prerequisite
    SOURCE_DIR "${SEQAN3_ROOT}"
    CMAKE_ARGS
        ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}
        "-DCMAKE_INSTALL_PREFIX=<BINARY_DIR>/usr"
    STEP_TARGETS configure install
    BUILD_BYPRODUCTS "<BINARY_DIR>/include"
)
set (SEQAN3_PACKAGE_ZIP_URL "${PROJECT_BINARY_DIR}/seqan3-${SEQAN3_VERSION}-${CMAKE_SYSTEM_NAME}.zip")
ExternalProject_Add_Step (
    seqan3_test_prerequisite
    package
    COMMAND ${CMAKE_CPACK_COMMAND} -G ZIP -B "${PROJECT_BINARY_DIR}"
    DEPENDEES configure install
    WORKING_DIRECTORY "<BINARY_DIR>"
    BYPRODUCTS
        ${SEQAN3_PACKAGE_ZIP_URL}
        "${SEQAN3_PACKAGE_ZIP_URL}.sha256"
)
ExternalProject_Get_property(seqan3_test_prerequisite BINARY_DIR)
set (SEQAN3_SYSTEM_PREFIX "${BINARY_DIR}/usr")
