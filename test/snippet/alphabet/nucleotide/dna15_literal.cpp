#include <seqan3/alphabet/nucleotide/dna15.hpp>

int main()
{
    using seqan3::operator""_dna15;
    
    seqan3::dna15_vector foo{"ACGTTA"_dna15};
    seqan3::dna15_vector bar = "ACGTTA"_dna15;
    auto bax = "ACGTTA"_dna15;
}
