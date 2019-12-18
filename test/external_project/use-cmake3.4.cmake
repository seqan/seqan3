# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.14)
project (seqan3_test_seqan3_config CXX)

# require Linux and x86_64
if (NOT (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_SIZEOF_VOID_P EQUAL 8))
    message (STATUS "cmake 3.4 test requires Linux and x86_64")
    return ()
endif ()

include (FetchContent)
FetchContent_Declare (
    cmake34
    URL "https://github.com/Kitware/CMake/releases/download/v3.4.3/cmake-3.4.3-Linux-x86_64.tar.gz"
    URL_HASH "SHA256=66b8d315c852908be9f79e1a18b8778714659fce4ddb2d041af8680a239202fc"
)
FetchContent_MakeAvailable(cmake34)

set (cmake34_command "${cmake34_SOURCE_DIR}/bin/cmake")

execute_process(COMMAND ${cmake34_command} --version
                RESULT_VARIABLE cmake34_result
                OUTPUT_VARIABLE cmake34_output)

if("${cmake34_result}" STREQUAL "0" AND cmake34_output MATCHES "cmake version 3.4.3")
    set (SEQAN3_EXTERNAL_PROJECT_CMAKE_COMMAND "${cmake34_command}")
    message (STATUS "Use cmake3.4 in tests [${cmake34_command}]")
else ()
    message (AUTHOR_WARNING "Couldn't execute cmake3.4 --version [${cmake34_command}]")
endif ()
