#include <seqan3/alphabet/nucleotide/dna16sam.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna16sam_vector sequence1{"ACGTTA"_dna16sam};
    seqan3::dna16sam_vector sequence2 = "ACGTTA"_dna16sam;
    auto sequence3 = "ACGTTA"_dna16sam;
}
