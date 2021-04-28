#include <seqan3/alphabet/quality/phred42.hpp>

int main()
{
    using namespace seqan3::literals;

    std::vector<seqan3::phred42> sequence1{"##!!##"_phred42};
    std::vector<seqan3::phred42> sequence2 = "##!!##"_phred42;
    auto sequence3 = "##!!##"_phred42;
}
