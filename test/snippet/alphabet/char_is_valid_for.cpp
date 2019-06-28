#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

int main()
{
    bool b = seqan3::char_is_valid_for<char>('A');
    // calls seqan3::custom::char_is_valid_for<char>('A'); always true

    bool c = seqan3::char_is_valid_for<seqan3::dna5>('A');
    // calls dna5::char_is_valid('A') member; == true

    // for some alphabets, characters that are not uniquely mappable are still valid:
    bool d = seqan3::char_is_valid_for<seqan3::dna5>('a');
    // lower case also true
}
