#include <iostream>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

int main()
{
    auto sigma_char = seqan3::alphabet_size<char>; // calls seqan3::custom::alphabet_size(char{})
    static_assert(std::same_as<decltype(sigma_char), uint16_t>);
    std::cout << sigma_char << '\n'; // 256

    auto sigma_dna5 = seqan3::alphabet_size<seqan3::dna5>; // returns dna5::alphabet_size
    static_assert(std::same_as<decltype(sigma_dna5), uint8_t>);
    std::cout << static_cast<uint16_t>(sigma_dna5) << '\n'; // 5
}
