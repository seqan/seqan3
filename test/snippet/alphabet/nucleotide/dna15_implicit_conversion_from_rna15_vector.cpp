// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

// generated from test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@_vector.cpp.in

//![main]
#include <vector>

#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna15_vector vector{'A'_rna15, 'C'_rna15, 'G'_rna15}; // (element-wise) implicit conversion

    // but this won't work:
    // seqan3::dna15_vector dna15_vector{"ACGT"_rna15};

    // as a workaround you can use:
    // side note: this would also work without the implicit conversion.
    seqan3::rna15_vector rna15_vector = "ACGT"_rna15;
    seqan3::dna15_vector dna15_vector{rna15_vector.begin(), rna15_vector.end()};
}
//![main]
