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

    // Default (first forward frame)
    auto v1 = vec | seqan3::views::translate_single;
    // == [T,Y,V,R]
    seqan3::debug_stream << v1[1] << '\n';

    // First forward frame
    auto v2 = vec | seqan3::views::translate_single(seqan3::translation_frames::forward_frame0);
    // == [T,Y,V,R]

    // First reverse frame
    auto v3 = vec | seqan3::views::translate_single(seqan3::translation_frames::reverse_frame0);
    // == [Y,V,R,T]

    // Second forward frame
    auto v4 = vec | seqan3::views::translate_single(seqan3::translation_frames::forward_frame1);
    // == [R,T,Y,V]

    // Second reverse frame
    auto v5 = vec | seqan3::views::translate_single(seqan3::translation_frames::reverse_frame1);
    // == [T,Y,V,R]

    // Third forward frame
    auto v6 = vec | seqan3::views::translate_single(seqan3::translation_frames::forward_frame2);
    // == [V,R,T]

    // Third reverse frame
    auto v7 = vec | seqan3::views::translate_single(seqan3::translation_frames::reverse_frame2);
    // == [R,T,Y]

    // function syntax
    auto v8 = seqan3::views::translate_single(vec, seqan3::translation_frames::forward_frame0);
    // == [T,Y,V,R]

    // combinability
    auto v9 =
        vec | seqan3::views::complement | seqan3::views::translate_single(seqan3::translation_frames::reverse_frame0);
    // == [M,H,A,C]

    // combinability with default parameter
    auto v10 = vec | seqan3::views::complement | seqan3::views::translate_single;
    // == [C,M,H,A]

    // combinability with default parameter
    auto v11 = vec | seqan3::views::complement | seqan3::views::translate_single();
    // == [C,M,H,A]
}
