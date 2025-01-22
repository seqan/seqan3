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

    std::vector<seqan3::dna4> genome{"AATAATAAC"_dna4};
    seqan3::fm_index index{genome}; // build the index

    auto cur = index.cursor(); // create a cursor
    // cur.cycle_back();                                            // cycle_back on begin() is undefined behaviour!
    cur.extend_right("AAC"_dna4);                                 // search the sequence "AAC"
    seqan3::debug_stream << cur.path_label(genome) << '\n';       // prints "AAC"
    seqan3::debug_stream << cur.last_rank() << '\n';              // prints 1
    seqan3::debug_stream << cur.query_length() << '\n';           // prints 3
    auto [left_bound, right_bound] = cur.suffix_array_interval(); // Get the half-open suffix array interval.
    seqan3::debug_stream << '[' << left_bound << ',' << right_bound << ")\n"; // prints "[7,8)"

    cur.cycle_back();                                       // search the sequence "AAT"
    seqan3::debug_stream << cur.path_label(genome) << '\n'; // prints "AAT"
    seqan3::debug_stream << cur.last_rank() << '\n';        // prints 3
    seqan3::debug_stream << cur.query_length() << '\n';     // prints 3
    auto interval = cur.suffix_array_interval();            // Get the half-open suffix array interval.
    seqan3::debug_stream << '[' << interval.begin_position << ',' << interval.end_position << ")\n"; // prints "[8,10)"

    cur.cycle_back();                                        // "cur" doesn't change because the rightmost char
                                                             // is already the largest dna4 char.
    seqan3::debug_stream << cur.path_label(genome) << '\n';  // prints "AAT"
    seqan3::debug_stream << cur.last_rank() << '\n';         // prints 3
    seqan3::debug_stream << cur.query_length() << '\n';      // prints 3
    auto && [lb, rb] = cur.suffix_array_interval();          // Get the half-open suffix array interval.
    seqan3::debug_stream << '[' << lb << ',' << rb << ")\n"; // prints "[8,10)"
}
