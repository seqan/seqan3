#include <seqan3/alphabet/quality/phred42.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::phred42 letter1{'!'_phred42};
    auto letter2 = '!'_phred42;
}
