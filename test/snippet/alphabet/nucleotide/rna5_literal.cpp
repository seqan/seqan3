#include <seqan3/alphabet/nucleotide/rna5.hpp>

int main()
{
    using seqan3::operator""_rna5;
    
    seqan3::rna5_vector foo{"ACGTTA"_rna5};
    seqan3::rna5_vector bar = "ACGTTA"_rna5;
    auto bax = "ACGTTA"_rna5;
}
