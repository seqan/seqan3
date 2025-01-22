// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/utility/range/to.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna4_vector vec{"ACGGTC"_dna4};
    auto vec_view2 = seqan3::views::complement(vec);

    // re-convert to container
    seqan3::dna4_vector complemented = vec_view2 | seqan3::ranges::to<seqan3::dna4_vector>();
    assert(complemented == "TGCCAG"_dna4);

    // also possible in one step
    seqan3::dna4_vector reversed = vec | std::views::reverse | seqan3::ranges::to<seqan3::dna4_vector>();
    assert(reversed == "CTGGCA"_dna4);
}
