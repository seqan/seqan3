#include <seqan3/alphabet/nucleotide/dna5.hpp>

int main()
{
    using seqan3::operator""_dna5;
    
    seqan3::dna5_vector foo{"ACGTTA"_dna5};
    seqan3::dna5_vector bar = "ACGTTA"_dna5;
    auto bax = "ACGTTA"_dna5;
}
