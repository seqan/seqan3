// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <algorithm>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/complement.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna4_vector vec{"ACGGTC"_dna4};

    // these are synonymous:
    auto vec_view1 = vec | seqan3::views::complement;
    auto vec_view2 = seqan3::views::complement(vec);
    std::cout << std::boolalpha << (std::ranges::equal(vec_view1, vec_view2)) << '\n'; // true
    // both views "behave" like a collection of the elements 'T', 'G', 'C', 'C', 'A', 'G'
    // but can be copied cheaply et cetera
}
