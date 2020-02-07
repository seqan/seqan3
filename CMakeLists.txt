# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

# This file provides functionality common to the different test modules used by
# SeqAn3. To build tests, run cmake on one of the sub-folders in this directory
# which contain a CMakeLists.txt.

cmake_minimum_required (VERSION 3.4)

find_path (SEQAN3_MODULE_PATH "seqan3-config.cmake" HINTS "${CMAKE_CURRENT_LIST_DIR}/build_system/")
list (APPEND CMAKE_MODULE_PATH "${SEQAN3_MODULE_PATH}")

include (seqan3-config-version)

if (CMAKE_VERSION VERSION_LESS 3.12)
    project (seqan3 LANGUAGES CXX VERSION "${PACKAGE_VERSION}")
else ()
    project (
        seqan3 LANGUAGES CXX VERSION "${PACKAGE_VERSION}"
        DESCRIPTION "SeqAn3 -- the modern C++ library for sequence analysis" # since cmake 3.9
        HOMEPAGE_URL "https://github.com/seqan/seqan3" # since cmake 3.12
    )
endif ()

find_package (SeqAn3 3.0 REQUIRED HINTS ${SEQAN3_MODULE_PATH})

option(INSTALL_SEQAN3 "Enable installation of seqan3. (Projects embedding seqan3 may want to turn this OFF.)" ON)

if (INSTALL_SEQAN3)
    include (seqan3-install)
    include (seqan3-package)
endif ()
