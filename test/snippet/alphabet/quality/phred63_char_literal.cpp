#include <seqan3/alphabet/quality/phred63.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::phred63 letter1{'!'_phred63};
    auto letter2 = '!'_phred63;
}
