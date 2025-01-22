// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/bloom_filter/bloom_filter.hpp>

int main()
{
    seqan3::bloom_filter bf{seqan3::bin_size{8192u}};

    // Insert the values `126`, `712` and `237` into the Bloom Filter.
    bf.emplace(126);
    bf.emplace(712);
    bf.emplace(237);
}
