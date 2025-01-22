// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <iostream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/to_rank.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    // The alphabet normally needs to be converted to char explicitly:
    std::cout << seqan3::to_char('C'_dna5) << '\n'; // prints 'C'

    // The debug_stream, on the other hand, does this automatically:
    seqan3::debug_stream << 'C'_dna5 << '\n'; // prints 'C'

    // The debug_stream can also print all types that model std::ranges::input_range:
    std::vector<seqan3::dna5> vec{"ACGT"_dna5};
    seqan3::debug_stream << vec << '\n'; // prints "ACGT"

    // ranges of non-alphabets are printed comma-separated:
    seqan3::debug_stream << (vec | seqan3::views::to_rank) << '\n'; // prints "[0,1,2,3]"
}
