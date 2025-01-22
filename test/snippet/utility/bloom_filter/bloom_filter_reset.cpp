// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/utility/bloom_filter/bloom_filter.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::bloom_filter bf{seqan3::bin_size{8192u}, seqan3::hash_function_count{2u}};

    auto const sequence1 = "ACTGACTGACTGATC"_dna4;
    auto const sequence2 = "GTGACTGACTGACTCG"_dna4;
    auto const sequence3 = "AAAAAAACGATCGACA"_dna4;
    auto kmers = seqan3::views::kmer_hash(seqan3::ungapped{5u});

    // Insert all 5-mers of sequence1
    for (auto && value : sequence1 | kmers)
        bf.emplace(value);

    // Insert all 5-mers of sequence3
    for (auto && value : sequence3 | kmers)
        bf.emplace(value);

    // Count all 5-mers of sequence2
    seqan3::debug_stream << bf.count(sequence2 | kmers) << '\n'; // 9

    // Reset the Bloom Filter
    bf.reset();

    // After reset, no 5-mers are found
    seqan3::debug_stream << bf.count(sequence2 | kmers) << '\n'; // 0
}
