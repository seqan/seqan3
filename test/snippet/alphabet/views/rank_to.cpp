// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/rank_to.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    std::vector<int> vec{0, 1, 3, 3, 3, 2, 0, 3, 0};
    seqan3::debug_stream << (vec | seqan3::views::rank_to<seqan3::dna4>) << '\n'; // ACTTTGATA
    seqan3::debug_stream << (vec | seqan3::views::rank_to<seqan3::dna5>) << '\n'; // ACNNNGANA
}
