# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.20...3.31)
project (seqan3_app CXX)

# --- helper scripts
include (../find-package-diagnostics.cmake)
file (SHA256 "${SEQAN3_PACKAGE_ZIP_URL}" SEQAN3_PACKAGE_ZIP_HASH)
message (STATUS "SEQAN3_PACKAGE_ZIP_URL: ${SEQAN3_PACKAGE_ZIP_URL}")
message (STATUS "SEQAN3_PACKAGE_ZIP_HASH: SHA256=${SEQAN3_PACKAGE_ZIP_HASH}")
# ---

# fetch seqan3 sources (requires >= cmake 3.14)
include (FetchContent)
FetchContent_Declare (
    seqan3
    URL "${SEQAN3_PACKAGE_ZIP_URL}" # change these values
    URL_HASH "SHA256=${SEQAN3_PACKAGE_ZIP_HASH}" # change these values
)
FetchContent_MakeAvailable (seqan3)

# add seqan3 to search path
list (APPEND CMAKE_PREFIX_PATH "${seqan3_SOURCE_DIR}")

# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)

# build app with seqan3
add_executable (hello_world ../src/hello_world.cpp)
target_link_libraries (hello_world seqan3::seqan3)
install (TARGETS hello_world)
