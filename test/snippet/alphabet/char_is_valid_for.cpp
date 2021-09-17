#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

int main()
{
    // calls seqan3::custom::char_is_valid_for<char>('A')
    std::cout << std::boolalpha << seqan3::char_is_valid_for<char>('A') << '\n'; // always 'true'

    // calls dna5::char_is_valid('A') member
    std::cout << std::boolalpha << seqan3::char_is_valid_for<seqan3::dna5>('A') << '\n'; // true

    // for some alphabets, characters that are not uniquely mappable are still valid:
    std::cout << std::boolalpha << seqan3::char_is_valid_for<seqan3::dna5>('a') << '\n'; // true
}
