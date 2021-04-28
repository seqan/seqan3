#include <seqan3/alphabet/aminoacid/aa10murphy.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::aa10murphy_vector sequence1{"ACGTTA"_aa10murphy};
    seqan3::aa10murphy_vector sequence2 = "ACGTTA"_aa10murphy;
    auto sequence3 = "ACGTTA"_aa10murphy;
}
