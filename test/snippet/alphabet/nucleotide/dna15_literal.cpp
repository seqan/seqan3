#include <seqan3/alphabet/nucleotide/dna15.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna15_vector sequence1{"ACGTTA"_dna15};
    seqan3::dna15_vector sequence2 = "ACGTTA"_dna15;
    auto sequence3 = "ACGTTA"_dna15;
}
