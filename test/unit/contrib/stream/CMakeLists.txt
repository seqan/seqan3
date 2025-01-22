# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

if (BZIP2_FOUND)
    seqan3_test (bz2_istream_test.cpp)
    seqan3_test (bz2_ostream_test.cpp)
endif ()

if (ZLIB_FOUND)
    seqan3_test (gz_istream_test.cpp)
    seqan3_test (gz_ostream_test.cpp)

    seqan3_test (bgzf_istream_test.cpp)
    seqan3_test (bgzf_ostream_test.cpp)
endif ()
