// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

// generated from test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@.cpp.in

//![main]
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna5 letter1 = 'C'_rna5; // implicitly converted
    seqan3::dna5 letter2{};
    letter2 = 'C'_rna5; // implicitly converted
}
//![main]
