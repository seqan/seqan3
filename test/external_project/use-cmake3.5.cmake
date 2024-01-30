# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.14)
project (seqan3_test_seqan3_config CXX)

# require Linux and x86_64
if (NOT (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_SIZEOF_VOID_P EQUAL 8))
    message (STATUS "cmake 3.5 test requires Linux and x86_64")
    return ()
endif ()

include (FetchContent)
FetchContent_Declare (
    cmake35
    URL "https://github.com/Kitware/CMake/releases/download/v3.5.2/cmake-3.5.2-Linux-x86_64.tar.gz"
    URL_HASH "SHA256=5f7aeaebe33521647625e0411467de71a2886743e4aa2c179e04c9e141c6c8cd")
FetchContent_MakeAvailable (cmake35)

set (cmake35_command "${cmake35_SOURCE_DIR}/bin/cmake")

execute_process (COMMAND ${cmake35_command} --version
                 RESULT_VARIABLE cmake35_result
                 OUTPUT_VARIABLE cmake35_output)

if ("${cmake35_result}" STREQUAL "0" AND cmake35_output MATCHES "cmake version 3.5.2")
    set (SEQAN3_EXTERNAL_PROJECT_CMAKE_COMMAND "${cmake35_command}")
    message (STATUS "Use cmake3.5 in tests [${cmake35_command}]")
else ()
    message (AUTHOR_WARNING "Couldn't execute cmake3.5 --version [${cmake35_command}]")
endif ()
