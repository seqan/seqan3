# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

option (HEADER_FILE_ABSOLUTE "")
option (HEADER_FILE_INCLUDE "")
option (HEADER_TARGET_SOURCE "")
option (HEADER_TEST_NAME_SAFE "")
option (HEADER_COMPONENT "")
option (HEADER_SUB_TEST "")

file (WRITE "${HEADER_TARGET_SOURCE}" "") # write empty file

# cmake-format: off

if (HEADER_SUB_TEST STREQUAL "no-self-include")
    # this test ensures that a header will not be included by itself later
    file (READ "${HEADER_FILE_ABSOLUTE}" header_content)

    string (REPLACE "#pragma once" "" header_content "${header_content}")

    file (APPEND "${HEADER_TARGET_SOURCE}"
        "// seqan3-header-test-no-self-include-start\n"
        "${header_content}\n"
        "// seqan3-header-test-no-self-include-end\n\n")
else ()
    # this test ensures that a header guard is in place
    file (APPEND "${HEADER_TARGET_SOURCE}" #
          "// seqan3-header-test-header-guard-start\n"
          "#include <${HEADER_FILE_INCLUDE}>\n" #
          "#include <${HEADER_FILE_INCLUDE}>\n"
          "// seqan3-header-test-header-guard-end\n\n")
endif ()

# these includes are required by some headers (note that they follow)
file (APPEND "${HEADER_TARGET_SOURCE}" #
      "// seqan3-header-test-dependencies-start\n"
      "#include <gtest/gtest.h>\n" #
      "#include <benchmark/benchmark.h>\n" #
      "TEST(${HEADER_TEST_NAME_SAFE}) {}\n"
      "// seqan3-header-test-dependencies-end\n\n")

# test that seqan3 headers include platform.hpp
if ("${HEADER_COMPONENT}" MATCHES "seqan3")

    # exclude seqan3/std/* and seqan3/contrib/* from platform test
    if (NOT HEADER_FILE_INCLUDE MATCHES "seqan3/(std|contrib)/")
        file (APPEND "${HEADER_TARGET_SOURCE}" #
              "// seqan3-header-test-platform-start\n"
              "#ifndef SEQAN3_DOXYGEN_ONLY\n" #
              "#error \"Your header '${HEADER_FILE_INCLUDE}' file is missing #include <seqan3/core/platform.hpp>\"\n" #
              "#endif\n"
              "// seqan3-header-test-platform-end\n\n")
    endif ()

    # seqan3/std/* must not include platform.hpp (and therefore any other seqan3 header)
    # See https://github.com/seqan/product_backlog/issues/135
    if (HEADER_FILE_INCLUDE MATCHES "seqan3/std/")
        file (APPEND "${HEADER_TARGET_SOURCE}" #
              "// seqan3-header-test-no-platform-start\n"
              "#ifdef SEQAN3_DOXYGEN_ONLY" #
              "#error \"The standard header '${HEADER_FILE_INCLUDE}' file MUST NOT include any other " #
              "seqan3 header (except for seqan3/contrib)\"\n" #
              "#endif\n"
              "// seqan3-header-test-no-platform-end\n\n")
    endif ()
endif ()

# cmake-format: on
