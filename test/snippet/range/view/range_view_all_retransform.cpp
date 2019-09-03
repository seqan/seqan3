#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/view/complement.hpp>
#include <seqan3/std/ranges>

int main()
{
    using seqan3::operator""_dna4;

    seqan3::dna4_vector vec{"ACGGTC"_dna4};
    auto vec_view2 = seqan3::view::complement(vec);

    // just re-assign to a container
    seqan3::dna4_vector complemented = vec_view2;
    assert(std::ranges::equal(complemented, "TGCCAG"_dna4));

    // or immediately create on container
    seqan3::dna4_vector reversed = vec | std::view::reverse;
    assert(std::ranges::equal(reversed, "CTGGCA"_dna4));
}
