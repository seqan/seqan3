#include <seqan3/alphabet/nucleotide/dna4.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna4 letter = 'A'_dna4; // identical to assign_char_to('A', letter);
    seqan3::dna4_vector sequence = "ACGT"_dna4; // identical to calling assign_char for each element
}
