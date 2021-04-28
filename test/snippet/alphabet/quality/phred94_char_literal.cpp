#include <seqan3/alphabet/quality/phred94.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::phred94 letter1{'!'_phred94};
    auto letter2 = '!'_phred94;
}
