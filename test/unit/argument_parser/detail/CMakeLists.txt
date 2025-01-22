# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

find_path (SEQAN3_TEST_LICENSE_DIR
           NAMES LICENSE.md
           HINTS "${SEQAN3_CLONE_DIR}")

add_definitions (-DSEQAN3_TEST_LICENSE_DIR="${SEQAN3_TEST_LICENSE_DIR}")

seqan3_test (format_help_test.cpp CYCLIC_DEPENDING_INCLUDES include-seqan3-argument_parser-detail-format_html.hpp
                                                            include-seqan3-argument_parser-detail-format_man.hpp)
seqan3_test (format_html_test.cpp CYCLIC_DEPENDING_INCLUDES include-seqan3-argument_parser-detail-format_help.hpp
                                                            include-seqan3-argument_parser-detail-format_man.hpp)
seqan3_test (format_man_test.cpp)
seqan3_test (version_check_debug_test.cpp)
seqan3_test (version_check_release_test.cpp)
seqan3_test (concept_is_container_option_test.cpp)
