#include <seqan3/alphabet/aminoacid/aa27.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::aa27_vector sequence1{"ACGTTA"_aa27};
    seqan3::aa27_vector sequence2 = "ACGTTA"_aa27;
    auto sequence3 = "ACGTTA"_aa27;
}
