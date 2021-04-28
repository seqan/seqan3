#include <seqan3/alphabet/nucleotide/dna4.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna4_vector sequence1{"ACGTTA"_dna4};
    seqan3::dna4_vector sequence2 = "ACGTTA"_dna4;
    auto sequence3 = "ACGTTA"_dna4;
}
