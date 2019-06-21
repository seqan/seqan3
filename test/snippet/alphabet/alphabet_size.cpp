#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

int main()
{
    auto r2 = seqan3::alphabet_size<char>;            // calls seqan3::custom::alphabet_size(char{}); r2 == 256
    auto r3 = seqan3::alphabet_size<seqan3::dna5>;    // returns dna5::alphabet_size; r3 == 5
}
