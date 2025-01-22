// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/bloom_filter/bloom_filter.hpp>

int main()
{
    // Construct a Bloom Filter that contains 8192 bits and 3 hash functions.
    seqan3::bloom_filter bf{seqan3::bin_size{8192u}, seqan3::hash_function_count{3}};

    // Construct a Bloom Filter that contains 256 KiBits and the default of 2
    // hash functions.
    seqan3::bloom_filter bf2{seqan3::bin_size{1ULL << 20}};
}
