#include <seqan3/alphabet/nucleotide/dna15.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna15 letter1{'A'_dna15};
    auto letter2 = 'A'_dna15;
}
