# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

add_subdirectories ()

seqan3_test (char_to_test.cpp)
seqan3_test (char_strictly_to_test.cpp)
seqan3_test (complement_test.cpp)
seqan3_test (rank_to_test.cpp)
seqan3_test (to_char_test.cpp)
seqan3_test (to_rank_test.cpp)
seqan3_test (translate_join_test.cpp)
seqan3_test (translate_test.cpp)
seqan3_test (trim_quality_test.cpp)
seqan3_test (validate_char_for_test.cpp)
