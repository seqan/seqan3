// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/bloom_filter/bloom_filter.hpp>

int main()
{
    // Construct an uncompressed Bloom Filter.
    seqan3::bloom_filter bf{seqan3::bin_size{128u}, seqan3::hash_function_count{3}};

    // Fill `bf` with data.

    // Construct an immutable, compressed Bloom Filter.
    seqan3::bloom_filter<seqan3::data_layout::compressed> bf2{bf};
}
