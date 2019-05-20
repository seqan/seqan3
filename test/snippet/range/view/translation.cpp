#include <iostream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/translation.hpp>
#include <seqan3/range/view/complement.hpp>

using namespace seqan3;

int main()
{
//! [dna5]
dna5_vector vec{"ACGTACGTACGTA"_dna5};

// Default (first forward frame)
auto v1 = vec | view::translate_single;                                                           // == [T,Y,V,R]
debug_stream << v1[1] << '\n';

// First forward frame
auto v2 = vec | view::translate_single(translation_frames::FWD_FRAME_0);                          // == [T,Y,V,R]

// First reverse frame
auto v3 = vec | view::translate_single(translation_frames::REV_FRAME_0);                          // == [Y,V,R,T]

// Second forward frame
auto v4 = vec | view::translate_single(translation_frames::FWD_FRAME_1);                          // == [R,T,Y,V]

// Second reverse frame
auto v5 = vec | view::translate_single(translation_frames::REV_FRAME_1);                          // == [T,Y,V,R]

// Third forward frame
auto v6 = vec | view::translate_single(translation_frames::FWD_FRAME_2);                            // == [V,R,T]

// Third reverse frame
auto v7 = vec | view::translate_single(translation_frames::REV_FRAME_2);                            // == [R,T,Y]

// function syntax
auto v8 = view::translate_single(vec, translation_frames::FWD_FRAME_0);                           // == [T,Y,V,R]

// combinability
auto v9 = vec | view::complement | view::translate_single(translation_frames::REV_FRAME_0);      // == [M,H,A,C]

// combinability with default parameter
auto v10 = vec | view::complement | view::translate_single;                                       // == [C,M,H,A]

// combinability with default parameter
auto v12 = vec | view::complement | view::translate_single();                                     // == [C,M,H,A]
//! [dna5]

{
//! [usage]
dna5_vector vec{"ACGTACGTACGTA"_dna5};

// default frame translation
auto v1 = vec | view::translate;                                  // == [[T,Y,V,R],[R,T,Y,V],[V,R,T],[Y,V,R,T],[T,Y,V,R],[R,T,Y]]

// single frame translation
auto v2 = vec | view::translate(translation_frames::FWD_FRAME_0);                                               // == [[T,Y,V,R]]

// reverse translation
auto v3 = vec | view::translate(translation_frames::FWD_REV_0);                                       // == [[T,Y,V,R],[Y,V,R,T]]

// forward frames translation
auto v4 = vec | view::translate(translation_frames::FWD);                                     // == [[T,Y,V,R],[R,T,Y,V],[V,R,T]]

// six frame translation
auto v5 = vec | view::translate(translation_frames::SIX_FRAME);   // == [[T,Y,V,R],[R,T,Y,V],[V,R,T],[Y,V,R,T],[T,Y,V,R],[R,T,Y]]
// function syntax
auto v6 = view::translate(vec, translation_frames::FWD_REV_0);                                        // == [[T,Y,V,R],[Y,V,R,T]]

// combinability
auto v7 = vec | view::complement | view::translate(translation_frames::FWD_REV_0);                    // == [[C,M,H,A],[M,H,A,C]]
//! [usage]
}
}
