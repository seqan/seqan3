// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;
    seqan3::dna4_vector s{"ACTTTGATAA"_dna4};
    using iterator = seqan3::dna4_vector::iterator;
    auto v1 = std::ranges::subrange<iterator, iterator>{std::ranges::begin(s) + 2, std::ranges::end(s)}
            | seqan3::views::to_char; // == "TTTGATAA"

    seqan3::debug_stream << v1 << '\n';
}
