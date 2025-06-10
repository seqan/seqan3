# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# We can only set one CMAKE_SYSTEM_PREFIX_PATH, i.e. it cannot be a list.
# Hence we need to reuse the SEQAN3_SYSTEM_PREFIX.
if (NOT DEFINED SEQAN3_SYSTEM_PREFIX)
    message (FATAL_ERROR "SEQAN3_SYSTEM_PREFIX is not defined. Did you include install-seqan3.cmake before this file?")
endif ()

# install and package sharg library
ExternalProject_Add (
    sharg_test_prerequisite
    PREFIX sharg_test_prerequisite
    URL "https://github.com/seqan/sharg-parser/archive/be113bcffe49c0d62cbd65a191820f05386aa8da.tar.gz"
    URL_HASH SHA256=d4a723f58865d0a737299d3c6bd85addc912add89b93e6aba43d4836ed4980d1
    CMAKE_ARGS ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS} #
               "-DCMAKE_INSTALL_PREFIX=${SEQAN3_SYSTEM_PREFIX}"
    STEP_TARGETS configure install
    BUILD_BYPRODUCTS "<BINARY_DIR>/include")

ExternalProject_Get_Property (sharg_test_prerequisite SOURCE_DIR)
set (SHARG_ROOT "${SOURCE_DIR}")
