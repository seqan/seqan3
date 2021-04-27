#include <seqan3/alphabet/quality/phred63.hpp>

int main()
{
    using namespace seqan3::literals;

    std::vector<seqan3::phred63> sequence1{"##!!##"_phred63};
    std::vector<seqan3::phred63> sequence2 = "##!!##"_phred63;
    auto sequence3 = "##!!##"_phred63;
}
