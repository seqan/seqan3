// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

//![start]
#include <ranges>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/debug_stream.hpp>

using seqan3::operator""_dna5;

//![start]
auto my_reverse_complement = std::views::reverse
                           | std::views::transform(
                                 [](auto const d)
                                 {
                                     return seqan3::complement(d);
                                 });

//![end]
int main()
{
    std::vector<seqan3::dna5> vec{"ACCAGATTA"_dna5};
    seqan3::debug_stream << vec << '\n'; // will print "ACCAGATTA"

    auto v = vec | my_reverse_complement;

    seqan3::debug_stream << v << '\n'; // prints "TAATCTGGT"
}
//![end]
