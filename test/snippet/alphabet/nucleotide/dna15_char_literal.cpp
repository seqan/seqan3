// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

// generated from test/snippet/alphabet/nucleotide/@target_alphabet@_char_literal.cpp.in

//![main]
#include <seqan3/alphabet/nucleotide/dna15.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna15 letter1{'A'_dna15};
    auto letter2 = 'A'_dna15;
}
//![main]
