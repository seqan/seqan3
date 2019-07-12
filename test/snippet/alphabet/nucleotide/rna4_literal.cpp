#include <seqan3/alphabet/nucleotide/rna4.hpp>

int main()
{
    using seqan3::operator""_rna4;
    
    seqan3::rna4_vector foo{"ACGTTA"_rna4};
    seqan3::rna4_vector bar = "ACGTTA"_rna4;
    auto bax = "ACGTTA"_rna4;
}
