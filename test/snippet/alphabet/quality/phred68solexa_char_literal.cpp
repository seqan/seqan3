#include <seqan3/alphabet/quality/phred68solexa.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::phred68solexa letter1{'!'_phred68solexa};
    auto letter2 = '!'_phred68solexa;
}
