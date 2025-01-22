# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# This file provides functionality common to the different test modules used by
# SeqAn3. To build tests, run cmake on one of the sub-folders in this directory
# which contain a CMakeLists.txt.

cmake_minimum_required (VERSION 3.20...3.31)

find_path (SEQAN3_MODULE_PATH "seqan3-config.cmake" HINTS "${CMAKE_CURRENT_LIST_DIR}/cmake/")
list (APPEND CMAKE_MODULE_PATH "${SEQAN3_MODULE_PATH}")

include (seqan3-config-version)

project (seqan3
         LANGUAGES CXX
         VERSION "${SEQAN3_PROJECT_VERSION}"
         DESCRIPTION "SeqAn3 -- the modern C++ library for sequence analysis"
         HOMEPAGE_URL "https://github.com/seqan/seqan3")

find_package (SeqAn3 3.0 REQUIRED HINTS ${SEQAN3_MODULE_PATH})

option (INSTALL_SEQAN3 "Enable installation of seqan3. (Projects embedding seqan3 may want to turn this OFF.)" ON)

if (INSTALL_SEQAN3)
    include (seqan3-install)
    include (seqan3-package)
endif ()
