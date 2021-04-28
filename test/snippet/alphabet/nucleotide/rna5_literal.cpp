#include <seqan3/alphabet/nucleotide/rna5.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::rna5_vector sequence1{"ACGTTA"_rna5};
    seqan3::rna5_vector sequence2 = "ACGTTA"_rna5;
    auto sequence3 = "ACGTTA"_rna5;
}
