#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/concept.hpp>

int main()
{
    auto r2 = seqan3::alphabet_size_v<char>;            // calls seqan3::adaptation::alphabet_size(char{}); r2 == 256
    auto r3 = seqan3::alphabet_size_v<seqan3::dna5>;    // returns dna5::value_size; r3 == 5
}
