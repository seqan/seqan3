// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/all.hpp>

int main()
{
    using namespace seqan3::literals;

    std::vector<std::vector<seqan3::dna4>> genomes{"ATCTGACGAAGGCTAGCTAGCTAAGGGA"_dna4,
                                                   "TAGCTGAAGCCATTGGCATCTGATCGGACT"_dna4,
                                                   "ACTGAGCTCGTC"_dna4,
                                                   "TGCATGCACCCATCGACTGACTG"_dna4,
                                                   "GTACGTACGTTACG"_dna4};

    seqan3::fm_index index{genomes}; // build the index

    auto cur = index.cursor();                                         // create a cursor
    cur.extend_right("CTGA"_dna4);                                     // search the pattern "CTGA"
    seqan3::debug_stream << "Number of hits: " << cur.count() << '\n'; // outputs: 5
    seqan3::debug_stream << "Positions in the genomes: ";
    for (auto && pos : cur.locate()) // outputs: (3,16) (2,1) (1,3) (0,2) (1,19)
        seqan3::debug_stream << pos << ' ';
    seqan3::debug_stream << '\n';
    return 0;
}
