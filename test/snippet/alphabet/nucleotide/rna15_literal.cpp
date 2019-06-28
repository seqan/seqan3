#include <seqan3/alphabet/nucleotide/rna15.hpp>

int main()
{
    using seqan3::operator""_rna15;
    
    seqan3::rna15_vector foo{"ACGTTA"_rna15};
    seqan3::rna15_vector bar = "ACGTTA"_rna15;
    auto bax = "ACGTTA"_rna15;
}
