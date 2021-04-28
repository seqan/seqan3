#include <seqan3/alphabet/quality/phred68solexa.hpp>

int main()
{
    using namespace seqan3::literals;

    std::vector<seqan3::phred68solexa> sequence1{"##!!##"_phred68solexa};
    std::vector<seqan3::phred68solexa> sequence2 = "##!!##"_phred68solexa;
    auto sequence3 = "##!!##"_phred68solexa;
}
