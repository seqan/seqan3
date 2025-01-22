// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

// generated from test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@_views.cpp.in

//![main]
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/utility/views/convert.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna4_vector vector = "ACG"_dna4;

    auto rna4_view = vector | seqan3::views::convert<seqan3::rna4>;

    for (auto && chr : rna4_view) // converts lazily on-the-fly
    {
        static_assert(std::same_as<decltype(chr), seqan3::rna4 &&>);
    }
}
//![main]
