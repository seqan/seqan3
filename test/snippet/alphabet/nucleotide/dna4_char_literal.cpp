#include <seqan3/alphabet/nucleotide/dna4.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna4 letter1{'A'_dna4};
    auto letter2 = 'A'_dna4;
}
