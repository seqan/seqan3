// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/structure/dot_bracket3.hpp>
#include <seqan3/alphabet/structure/structured_rna.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;
    using seqan3::get;

    seqan3::structured_rna<seqan3::rna4, seqan3::dot_bracket3> letter{'G'_rna4, '('_db3};

    seqan3::debug_stream << seqan3::to_rank(letter) << ' '           // 7
                         << seqan3::to_rank(get<0>(letter)) << ' '   // 2
                         << seqan3::to_rank(get<1>(letter)) << '\n'; // 1

    seqan3::debug_stream << seqan3::to_char(letter) << ' '           // G
                         << seqan3::to_char(get<0>(letter)) << ' '   // G
                         << seqan3::to_char(get<1>(letter)) << '\n'; // (

    // Modify:
    get<0>(letter) = 'U'_rna4;
    seqan3::debug_stream << seqan3::to_char(letter) << '\n'; // U
}
