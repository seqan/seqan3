#include <seqan3/alphabet/nucleotide/dna3bs.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna3bs_vector sequence1{"ACGTTA"_dna3bs};
    seqan3::dna3bs_vector sequence2 = "ACGTTA"_dna3bs;
    auto sequence3 = "ACGTTA"_dna3bs;
}
