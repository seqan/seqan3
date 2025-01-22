// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/alphabet/views/translate.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna5_vector vec{"ACGTACGTACGTA"_dna5};

    // default frame translation
    auto v1 = vec | seqan3::views::translate;
    seqan3::debug_stream << v1 << '\n'; // [TYVR,RTYV,VRT,YVRT,TYVR,RTY]

    // single frame translation
    auto v2 = vec | seqan3::views::translate(seqan3::translation_frames::forward_frame0);
    seqan3::debug_stream << v2 << '\n'; // [TYVR]

    // reverse translation
    auto v3 = vec | seqan3::views::translate(seqan3::translation_frames::forward_reverse0);
    seqan3::debug_stream << v3 << '\n'; // [TYVR,YVRT]

    // forward frames translation
    auto v4 = vec | seqan3::views::translate(seqan3::translation_frames::forward_frames);
    seqan3::debug_stream << v4 << '\n'; // [TYVR,RTYV,VRT]

    // six frame translation
    auto v5 = vec | seqan3::views::translate(seqan3::translation_frames::six_frames);
    seqan3::debug_stream << v5 << '\n'; // [TYVR,RTYV,VRT,YVRT,TYVR,RTY]

    // function syntax
    auto v6 = seqan3::views::translate(vec, seqan3::translation_frames::forward_reverse0);
    seqan3::debug_stream << v6 << '\n'; // [TYVR,YVRT]

    // combinability
    auto v7 = vec | seqan3::views::complement | seqan3::views::translate(seqan3::translation_frames::forward_reverse0);
    seqan3::debug_stream << v7 << '\n'; // [CMHA,MHAC]
}
