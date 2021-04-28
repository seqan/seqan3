#include <seqan3/alphabet/aminoacid/aa10li.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::aa10li_vector sequence1{"ACGTTA"_aa10li};
    seqan3::aa10li_vector sequence2 = "ACGTTA"_aa10li;
    auto sequence3 = "ACGTTA"_aa10li;
}
