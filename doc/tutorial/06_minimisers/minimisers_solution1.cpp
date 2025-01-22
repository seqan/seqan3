// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

using namespace seqan3::literals;

int main()
{
    std::vector<seqan3::dna4> text{"CCACGTCGACGGTT"_dna4};

    // Here a consecutive shape with size 4 (so the k-mer size is 4) and a window size of 8 is used.
    auto minimisers = text | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{4}}, seqan3::window_size{8});
    // results in: [10322096095657499224, 10322096095657499142, 10322096095657499224]
    // representing the k-mers [GACG, TCGA, GACG]
    seqan3::debug_stream << minimisers << '\n';
}
