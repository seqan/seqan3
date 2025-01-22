// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

int main()
{
    char c = '!';
    seqan3::assign_rank_to(66, c); // calls seqan3::custom::assign_rank_to(66, c); == 'B'

    seqan3::dna5 d{};
    seqan3::assign_rank_to(2, d); // calls .assign_rank(2) member; == 'G'_dna5

    // also works for temporaries:
    seqan3::dna5 d2 = seqan3::assign_rank_to(2, seqan3::dna5{});

    // too-large ranks are undefined behaviour:
    // seqan3::dna5 d3 = seqan3::assign_rank_to(50, seqan3::dna5{});
}
