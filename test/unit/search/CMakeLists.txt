# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

add_subdirectories ()

seqan3_test (sdsl_index_test.cpp)

seqan3_test (search_collection_test.cpp)
seqan3_test (search_configuration_test.cpp)
seqan3_test (search_scheme_algorithm_test.cpp)
seqan3_test (search_scheme_test.cpp)
seqan3_test (search_test.cpp)
