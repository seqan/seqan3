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
    using seqan3::get;

    seqan3::qualified<seqan3::dna4, seqan3::phred42> letter{'A'_dna4, '('_phred42};

    seqan3::debug_stream << seqan3::to_rank(letter) << ' '           // 7
                         << seqan3::to_rank(get<0>(letter)) << ' '   // 0
                         << seqan3::to_rank(get<1>(letter)) << '\n'; // 7

    seqan3::debug_stream << seqan3::to_char(letter) << ' '           // A
                         << seqan3::to_char(get<0>(letter)) << ' '   // A
                         << seqan3::to_char(get<1>(letter)) << '\n'; // (

    seqan3::debug_stream << seqan3::to_phred(letter) << ' '           // 7
                         << seqan3::to_phred(get<1>(letter)) << '\n'; // 7

    // Modify:
    get<0>(letter) = 'G'_dna4;
    seqan3::debug_stream << seqan3::to_char(letter) << '\n'; // G
}
