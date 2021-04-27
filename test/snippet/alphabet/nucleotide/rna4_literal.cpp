#include <seqan3/alphabet/nucleotide/rna4.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::rna4_vector sequence1{"ACGTTA"_rna4};
    seqan3::rna4_vector sequence2 = "ACGTTA"_rna4;
    auto sequence3 = "ACGTTA"_rna4;
}
