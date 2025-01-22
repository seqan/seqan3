# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.20...3.31)
project (seqan3_test_analyze CXX)

include (../seqan3-test.cmake)

enable_testing ()

option (SEQAN3_GENERATE_SNIPPETS "" OFF)

add_subdirectory (../performance performance)
add_subdirectory (../snippet snippet)
add_subdirectory (../unit unit)
