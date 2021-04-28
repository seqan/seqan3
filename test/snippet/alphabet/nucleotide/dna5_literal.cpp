#include <seqan3/alphabet/nucleotide/dna5.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna5_vector sequence1{"ACGTTA"_dna5};
    seqan3::dna5_vector sequence2 = "ACGTTA"_dna5;
    auto sequence3 = "ACGTTA"_dna5;
}
