cmake_minimum_required (VERSION 3.7)

option (HEADER_FILE_ABSOLUTE "")
option (HEADER_FILE_INCLUDE "")
option (HEADER_TARGET_SOURCE "")
option (HEADER_TEST_NAME_SAFE "")
option (HEADER_COMPONENT "")
option (HEADER_SUB_TEST "")

file (WRITE "${HEADER_TARGET_SOURCE}" "") # write empty file

if (HEADER_SUB_TEST STREQUAL "no-self-include")
    # this test ensures that a header will not be included by itself later
    file (READ "${HEADER_FILE_ABSOLUTE}" header_content)

    string (REPLACE "#pragma once" "" header_content "${header_content}")

    file (APPEND "${HEADER_TARGET_SOURCE}" "${header_content}")
else ()
    # this test ensures that a header guard is in place
    file (APPEND "${HEADER_TARGET_SOURCE}" "
#include <${HEADER_FILE_INCLUDE}>
#include <${HEADER_FILE_INCLUDE}>")
endif()

# these includes are required by some headers (note that they follow)
file (APPEND "${HEADER_TARGET_SOURCE}" "
#include <gtest/gtest.h>
#include <benchmark/benchmark.h>
TEST(${HEADER_TEST_NAME_SAFE}) {}")

# test that seqan3 headers include platform.hpp
if ("${HEADER_COMPONENT}" MATCHES "seqan3")

    # exclude seqan3/std/* and seqan3/contrib/* from platform test
    if (NOT HEADER_FILE_INCLUDE MATCHES "seqan3/(std|contrib)/")
        file (APPEND "${HEADER_TARGET_SOURCE}" "
#ifndef SEQAN3_DOXYGEN_ONLY
#error \"Your header '${HEADER_FILE_INCLUDE}' file is missing #include <seqan3/core/platform.hpp>\"
#endif")
    endif ()

    # seqan3/std/* must not include platform.hpp (and therefore any other seqan3 header)
    # See https://github.com/seqan/product_backlog/issues/135
    if (HEADER_FILE_INCLUDE MATCHES "seqan3/std/")
        file (APPEND "${HEADER_TARGET_SOURCE}" "
#ifdef SEQAN3_DOXYGEN_ONLY
#error \"The standard header '${HEADER_FILE_INCLUDE}' file MUST NOT include any other seqan3 header (except for seqan3/contrib)\"
#endif")
    endif ()

    # test whether seqan3 has the visibility bug on lower gcc versions
    # https://github.com/seqan/seqan3/issues/1317
    file (APPEND "${HEADER_TARGET_SOURCE}" "
#include <seqan3/core/platform.hpp>
class A{ int i{5}; };

template <typename t>
SEQAN3_CONCEPT private_bug = requires(t a){a.i;};

static_assert(!private_bug<A>, \"See https://github.com/seqan/seqan3/issues/1317\");")
endif ()
