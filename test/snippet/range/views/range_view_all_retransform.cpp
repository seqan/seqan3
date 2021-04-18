#include <seqan3/std/ranges>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/utility/views/to.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna4_vector vec{"ACGGTC"_dna4};
    auto vec_view2 = seqan3::views::complement(vec);

    // re-convert to container
    seqan3::dna4_vector complemented = vec_view2 | seqan3::views::to<seqan3::dna4_vector>;
    assert(complemented == "TGCCAG"_dna4);

    // also possible in one step
    seqan3::dna4_vector reversed = vec | std::views::reverse | seqan3::views::to<seqan3::dna4_vector>;
    assert(reversed == "CTGGCA"_dna4);
}
