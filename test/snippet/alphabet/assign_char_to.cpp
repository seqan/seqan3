// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

int main()
{
    char c = '!';
    seqan3::assign_char_to('?', c); // calls seqan3::custom::assign_char_to('A', c)

    seqan3::dna5 d{};
    seqan3::assign_char_to('A', d); // calls .assign_char('A') member

    // also works for temporaries:
    seqan3::dna5 d2 = seqan3::assign_char_to('A', seqan3::dna5{});

    // invalid/unknown characters are converted:
    seqan3::dna5 d3 = seqan3::assign_char_to('!', seqan3::dna5{}); // == 'N'_dna5
}
