#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/alphabet/views/translate.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna5_vector vec{"ACGTACGTACGTA"_dna5};

    // default frame translation
    auto v1 = vec | seqan3::views::translate;
    // == [[T,Y,V,R],[R,T,Y,V],[V,R,T],[Y,V,R,T],[T,Y,V,R],[R,T,Y]]

    // single frame translation
    auto v2 = vec | seqan3::views::translate(seqan3::translation_frames::forward_frame0);
    // == [[T,Y,V,R]]

    // reverse translation
    auto v3 = vec | seqan3::views::translate(seqan3::translation_frames::forward_reverse0);
    // == [[T,Y,V,R],[Y,V,R,T]]

    // forward frames translation
    auto v4 = vec | seqan3::views::translate(seqan3::translation_frames::forward_frames);
    // == [[T,Y,V,R],[R,T,Y,V],[V,R,T]]

    // six frame translation
    auto v5 = vec | seqan3::views::translate(seqan3::translation_frames::six_frames);
    // == [[T,Y,V,R],[R,T,Y,V],[V,R,T],[Y,V,R,T],[T,Y,V,R],[R,T,Y]]

    // function syntax
    auto v6 = seqan3::views::translate(vec, seqan3::translation_frames::forward_reverse0);
    // == [[T,Y,V,R],[Y,V,R,T]]

    // combinability
    auto v7 = vec | seqan3::views::complement | seqan3::views::translate(seqan3::translation_frames::forward_reverse0);
    // == [[C,M,H,A],[M,H,A,C]]
}
