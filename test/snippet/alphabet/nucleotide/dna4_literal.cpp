#include <seqan3/alphabet/nucleotide/dna4.hpp>

int main()
{
    using seqan3::operator""_dna4;
    
    seqan3::dna4_vector foo{"ACGTTA"_dna4};
    seqan3::dna4_vector bar = "ACGTTA"_dna4;
    auto bax = "ACGTTA"_dna4;
}
