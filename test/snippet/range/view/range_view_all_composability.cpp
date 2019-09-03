#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/view/complement.hpp>
#include <seqan3/std/ranges>

int main()
{
    using seqan3::operator""_dna4;

    seqan3::dna4_vector vec{"ACGGTC"_dna4};
    // views can be composed iteratively
    auto vec_view3 = vec | std::view::reverse;
    auto vec_view4 = vec_view3 | seqan3::view::complement;

    // or in one line similar to the unix command line
    auto vec_view5 = vec | seqan3::view::complement | std::view::reverse;

    // vec_view4 and vec_view5 are the reverse complement of "ACGGTC": "GACCGT"
}
