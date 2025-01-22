// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/structure/dssp9.hpp>
#include <seqan3/alphabet/structure/structured_aa.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;
    using seqan3::get;

    seqan3::structured_aa<seqan3::aa27, seqan3::dssp9> letter{'W'_aa27, 'B'_dssp9};

    seqan3::debug_stream << seqan3::to_rank(letter) << ' '           // 199
                         << seqan3::to_rank(get<0>(letter)) << ' '   // 22
                         << seqan3::to_rank(get<1>(letter)) << '\n'; // 1

    seqan3::debug_stream << seqan3::to_char(letter) << ' '           // W
                         << seqan3::to_char(get<0>(letter)) << ' '   // W
                         << seqan3::to_char(get<1>(letter)) << '\n'; // B

    // Modify:
    get<0>(letter) = 'V'_aa27;
    seqan3::debug_stream << seqan3::to_char(letter) << '\n'; // V
}
