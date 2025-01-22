// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

using namespace seqan3::literals;

int main()
{
    std::vector<seqan3::dna4> text{"ACGTAGC"_dna4};

    auto hashes = text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{3}});
    seqan3::debug_stream << hashes << '\n'; // [6,27,44,50,9]

    seqan3::debug_stream << (text | seqan3::views::kmer_hash(seqan3::ungapped{3})) << '\n'; // [6,27,44,50,9]

    seqan3::debug_stream << (text | seqan3::views::kmer_hash(0b101_shape)) << '\n'; // [2,7,8,14,1]

    // Attention: the Shape is defined from right to left!
    // The mask 0b11111101 applied to "AGAAAATA" ("A.AAAATA") will yield
    // the same hash value as mask 0b1111111 applied to "AAAAATA".
    {
        auto text1 = "AGAAAATA"_dna4;
        auto text2 = "AAAAATA"_dna4;
        seqan3::debug_stream << (text1 | seqan3::views::kmer_hash(0b1111'1101_shape)) << '\n'; // [12]
        seqan3::debug_stream << (text2 | seqan3::views::kmer_hash(0b111'1111_shape)) << '\n';  // [12]
    }
}
