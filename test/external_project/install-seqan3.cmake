# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# install and package seqan3 library
ExternalProject_Add (
    seqan3_test_prerequisite
    PREFIX seqan3_test_prerequisite
    SOURCE_DIR "${SEQAN3_ROOT}"
    CMAKE_ARGS ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS} #
               "-DCMAKE_INSTALL_PREFIX=<BINARY_DIR>/usr"
    STEP_TARGETS configure install
    BUILD_BYPRODUCTS "<BINARY_DIR>/include")

set (SEQAN3_PACKAGE_ZIP_URL "${PROJECT_BINARY_DIR}/seqan3-${SEQAN3_VERSION}-${CMAKE_SYSTEM_NAME}.zip")
ExternalProject_Add_Step (
    seqan3_test_prerequisite package
    COMMAND ${CMAKE_CPACK_COMMAND} -G ZIP -B "${PROJECT_BINARY_DIR}"
    DEPENDEES configure install
    WORKING_DIRECTORY "<BINARY_DIR>"
    BYPRODUCTS ${SEQAN3_PACKAGE_ZIP_URL} #
               "${SEQAN3_PACKAGE_ZIP_URL}.sha256")

ExternalProject_Get_Property (seqan3_test_prerequisite BINARY_DIR)
set (SEQAN3_SYSTEM_PREFIX "${BINARY_DIR}/usr")
