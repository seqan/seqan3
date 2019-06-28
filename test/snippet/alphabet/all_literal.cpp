#include <seqan3/alphabet/nucleotide/dna4.hpp>

int main()
{
    using seqan3::operator""_dna4;
    seqan3::dna4        my_letter = 'A'_dna4;           // identical to assign_char_to('A', my_letter);
    seqan3::dna4_vector my_seq    = "ACGT"_dna4;        // identical to calling assign_char for each element
}
