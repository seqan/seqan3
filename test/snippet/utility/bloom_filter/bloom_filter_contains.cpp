// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/bloom_filter/bloom_filter.hpp>

int main()
{
    seqan3::bloom_filter bf{seqan3::bin_size{1024u}};
    bf.emplace(126);
    bf.emplace(712);
    bf.emplace(237);

    // Query the Bloom Filter.
    // A return of `false` guarantees the query not being present in the Bloom Filter.
    // A return of `true` indicates the (probable) presence of the query in the Bloom Filter.
    // Note that there may be false positive results!
    bool result = bf.contains(712);
    seqan3::debug_stream << result << '\n'; // prints 1
}
