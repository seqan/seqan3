// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna5_vector foo{"ACGTA"_dna5};

    // pipe notation
    auto v = foo | seqan3::views::complement;
    seqan3::debug_stream << v << '\n'; // TGCAT

    // function notation
    auto v2(seqan3::views::complement(foo));
    seqan3::debug_stream << v2 << '\n'; // TGCAT

    // generate the reverse complement:
    auto v3 = foo | seqan3::views::complement | std::views::reverse;
    seqan3::debug_stream << v3 << '\n'; // TACGT
}
