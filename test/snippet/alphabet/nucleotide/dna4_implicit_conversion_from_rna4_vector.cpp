// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

// generated from test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@_vector.cpp.in

//![main]
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna4_vector vector{'A'_rna4, 'C'_rna4, 'G'_rna4}; // (element-wise) implicit conversion

    // but this won't work:
    // seqan3::dna4_vector dna4_vector{"ACGT"_rna4};

    // as a workaround you can use:
    // side note: this would also work without the implicit conversion.
    seqan3::rna4_vector rna4_vector = "ACGT"_rna4;
    seqan3::dna4_vector dna4_vector{rna4_vector.begin(), rna4_vector.end()};
}
//![main]
