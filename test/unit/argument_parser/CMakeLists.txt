# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

add_subdirectories ()
seqan3_test (format_parse_test.cpp)
seqan3_test (format_parse_validators_test.cpp)
seqan3_test (argument_parser_design_error_test.cpp)
