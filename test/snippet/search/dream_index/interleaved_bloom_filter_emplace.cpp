// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

int main()
{
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{12u}, seqan3::bin_size{8192u}};

    // Insert the values `126`, `712` and `237` into bins `0`, `3` and `9` of the Interleaved Bloom Filter.
    ibf.emplace(126, seqan3::bin_index{0u});
    ibf.emplace(712, seqan3::bin_index{3u});
    ibf.emplace(237, seqan3::bin_index{9u});
}
