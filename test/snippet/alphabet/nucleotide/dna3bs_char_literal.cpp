#include <seqan3/alphabet/nucleotide/dna3bs.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna3bs letter1{'A'_dna3bs};
    auto letter2 = 'A'_dna3bs;
}
