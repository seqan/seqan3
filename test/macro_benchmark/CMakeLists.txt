# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.10...3.30)
project (seqan3_test_unit CXX)

include (../seqan3-test.cmake)

include (seqan3_macro_benchmark)

seqan3_require_ccache ()
seqan3_require_benchmark ()

add_subdirectories ()
