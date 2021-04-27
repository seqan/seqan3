#include <seqan3/alphabet/nucleotide/rna15.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::rna15_vector sequence1{"ACGTTA"_rna15};
    seqan3::rna15_vector sequence2 = "ACGTTA"_rna15;
    auto sequence3 = "ACGTTA"_rna15;
}
