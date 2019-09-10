//! [example]
#include <iostream>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/char_to.hpp>
#include <seqan3/range/views/translate_join.hpp>

using seqan3::operator""_dna4;

int main()
{
    // Input range needs to be two-dimensional
    std::vector<std::vector<seqan3::dna4> > vec{"ACGTACGTACGTA"_dna4, "TCGAGAGCTTTAGC"_dna4};

    // Translation with default parameters
    auto v1 = vec | seqan3::views::translate_join;
    seqan3::debug_stream << v1 << "\n"; // [TYVR,RTYV,VRT,YVRT,TYVR,RTY,SRAL,REL*,ESFS,AKAL,LKLS,*SSR]

    // Access the third forward frame (index_frame 2) of the second input sequence (index_seq 1)
    // Required frames per sequence s = 6
    // n = (index_seq * s) + j
    //   = 1 * 6 + 2
    //   = 8

    auto third_frame_second_seq = v1[1 * 6 + 2];
    seqan3::debug_stream << third_frame_second_seq << "\n"; // ESFS

    // Translation with custom translation frame
    auto v2 = vec | seqan3::views::translate_join(seqan3::translation_frames::FWD_FRAME_0);
    seqan3::debug_stream << v2 << "\n"; // [TYVR,SRAL]

    return 0;
}
//! [example]
