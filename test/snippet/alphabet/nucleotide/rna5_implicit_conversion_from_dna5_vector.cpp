// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

// generated from test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@_vector.cpp.in

//![main]
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::rna5_vector vector{'A'_dna5, 'C'_dna5, 'G'_dna5}; // (element-wise) implicit conversion

    // but this won't work:
    // seqan3::rna5_vector rna5_vector{"ACGT"_dna5};

    // as a workaround you can use:
    // side note: this would also work without the implicit conversion.
    seqan3::dna5_vector dna5_vector = "ACGT"_dna5;
    seqan3::rna5_vector rna5_vector{dna5_vector.begin(), dna5_vector.end()};
}
//![main]
