# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.5)
project (seqan3_app CXX)

# --- helper scripts
include (../find-package-diagnostics.cmake)
# ---

# add seqan3 to search path
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/submodules/seqan3/build_system") # change this value

# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)

# build app with seqan3
add_executable (hello_world ../src/hello_world.cpp)
target_link_libraries (hello_world seqan3::seqan3)
if (CMAKE_VERSION VERSION_LESS 3.14)
    install (TARGETS hello_world RUNTIME DESTINATION bin)
else ()
    install (TARGETS hello_world) # RUNTIME DESTINATION not needed anymore since cmake 3.14
endif ()
