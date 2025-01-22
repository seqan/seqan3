// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <iostream>
#include <ranges> // include all of the standard library's views
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser.hpp>

using namespace seqan3::literals;

int main()
{
    std::vector<seqan3::dna4> text{"CCACGTCGACGGTT"_dna4};
    {
        //![minimiser]
        //![def]
        auto minimisers = text | seqan3::views::kmer_hash(seqan3::ungapped{4}) | seqan3::views::minimiser(5);
        //![minimiser]
    }

    {
        //![minimiser_seed]
        //![def]
        uint64_t const seed = 0x8F'3F'73'B5'CF'1C'9A'DE;
        auto minimisers = text | seqan3::views::kmer_hash(seqan3::ungapped{4})
                        | std::views::transform(
                              [](uint64_t i)
                              {
                                  return i ^ seed;
                              })
                        | seqan3::views::minimiser(5);
        //![minimiser_seed]
    }
}
