#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/complement.hpp>
#include <seqan3/range/views/translate.hpp>

int main()
{
    using seqan3::operator""_dna5;

    seqan3::dna5_vector vec{"ACGTACGTACGTA"_dna5};

    // Default (first forward frame)
    auto v1 = vec | seqan3::views::translate_single;
    // == [T,Y,V,R]
    seqan3::debug_stream << v1[1] << '\n';

    // First forward frame
    auto v2 = vec | seqan3::views::translate_single(seqan3::translation_frames::FWD_FRAME_0);
    // == [T,Y,V,R]

    // First reverse frame
    auto v3 = vec | seqan3::views::translate_single(seqan3::translation_frames::REV_FRAME_0);
    // == [Y,V,R,T]

    // Second forward frame
    auto v4 = vec | seqan3::views::translate_single(seqan3::translation_frames::FWD_FRAME_1);
    // == [R,T,Y,V]

    // Second reverse frame
    auto v5 = vec | seqan3::views::translate_single(seqan3::translation_frames::REV_FRAME_1);
    // == [T,Y,V,R]

    // Third forward frame
    auto v6 = vec | seqan3::views::translate_single(seqan3::translation_frames::FWD_FRAME_2);
    // == [V,R,T]

    // Third reverse frame
    auto v7 = vec | seqan3::views::translate_single(seqan3::translation_frames::REV_FRAME_2);
    // == [R,T,Y]

    // function syntax
    auto v8 = seqan3::views::translate_single(vec, seqan3::translation_frames::FWD_FRAME_0);
    // == [T,Y,V,R]

    // combinability
    auto v9 = vec | seqan3::views::complement | seqan3::views::translate_single(seqan3::translation_frames::REV_FRAME_0);
    // == [M,H,A,C]

    // combinability with default parameter
    auto v10 = vec | seqan3::views::complement | seqan3::views::translate_single;
    // == [C,M,H,A]

    // combinability with default parameter
    auto v11 = vec | seqan3::views::complement | seqan3::views::translate_single();
    // == [C,M,H,A]
}
