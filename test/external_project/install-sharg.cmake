# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.14)

# https://cmake.org/cmake/help/latest/policy/CMP0135.html
if (POLICY CMP0135)
    cmake_policy (SET CMP0135 NEW)
endif ()

# We can only set one CMAKE_SYSTEM_PREFIX_PATH, i.e. it cannot be a list.
# Hence we need to reuse the SEQAN3_SYSTEM_PREFIX.
if (NOT DEFINED SEQAN3_SYSTEM_PREFIX)
    message (FATAL_ERROR "SEQAN3_SYSTEM_PREFIX is not defined. Did you include install-seqan3.cmake before this file?")
endif ()

# install and package sharg library
ExternalProject_Add (
    sharg_test_prerequisite
    PREFIX sharg_test_prerequisite
    URL "https://github.com/seqan/sharg-parser/releases/download/1.1.1/sharg-1.1.1-Source.tar.xz"
    URL_HASH SHA256=7330f06501718e7871e55e5fd70d0e41472cc8b34bd0e3519f8c5547510c671c
    CMAKE_ARGS ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS} #
               "-DCMAKE_INSTALL_PREFIX=${SEQAN3_SYSTEM_PREFIX}"
    STEP_TARGETS configure install
    BUILD_BYPRODUCTS "<BINARY_DIR>/include")

ExternalProject_Get_Property (sharg_test_prerequisite SOURCE_DIR)
set (SHARG_ROOT "${SOURCE_DIR}")
