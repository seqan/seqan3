// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::qualified<seqan3::dna4, seqan3::phred42> letter1{'C'_dna4};    // creates {'C'_dna4, '!'_phred42}
    seqan3::qualified<seqan3::dna4, seqan3::phred42> letter2{'"'_phred42}; // creates {'A'_dna4, '"'_phred42}

    if (letter1 == letter2)
        seqan3::debug_stream << "yeah\n"; // yeah
}
