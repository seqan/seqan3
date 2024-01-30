// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

// generated from test/snippet/alphabet/nucleotide/@target_alphabet@_literal.cpp.in

//![main]
#include <seqan3/alphabet/nucleotide/dna3bs.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna3bs_vector sequence1{"ACGTTA"_dna3bs};
    seqan3::dna3bs_vector sequence2 = "ACGTTA"_dna3bs;
    auto sequence3 = "ACGTTA"_dna3bs;
}
//![main]
