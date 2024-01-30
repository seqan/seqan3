// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::qualified<seqan3::dna4, seqan3::phred42> letter1{'T'_dna4, '"'_phred42};

    letter1 = 'C'_dna4;    // yields {'C'_dna4, '"'_phred42}
    letter1 = '#'_phred42; // yields {'C'_dna4, '#'_phred42}
}
