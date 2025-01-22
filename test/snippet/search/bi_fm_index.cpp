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
    std::vector<seqan3::dna4> genome{"ATCGATCGAAGGCTAGCTAGCTAAGGGA"_dna4};
    seqan3::bi_fm_index index{genome}; // build the index

    auto cur = index.cursor();                                         // create a cursor
    cur.extend_right("GG"_dna4);                                       // search the pattern "GG"
    cur.extend_left("AA"_dna4);                                        // search the pattern "AAGG"
    seqan3::debug_stream << "Number of hits: " << cur.count() << '\n'; // outputs: 2
    seqan3::debug_stream << "Positions in the genome: ";
    for (auto && pos : cur.locate()) // outputs: (0, 8), (0, 22)
        seqan3::debug_stream << pos << ' ';
    seqan3::debug_stream << '\n';
    return 0;
}
