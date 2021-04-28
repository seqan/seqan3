#include <seqan3/alphabet/cigar/cigar.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::cigar::operation letter1{'M'_cigar_operation};
    auto letter2 = 'M'_cigar_operation;
}
