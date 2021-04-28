#include <seqan3/alphabet/nucleotide/dna16sam.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna16sam letter1{'A'_dna16sam};
    auto letter2 = 'A'_dna16sam;
}
