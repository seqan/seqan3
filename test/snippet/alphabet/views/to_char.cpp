#include <iostream>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

int main()
{
    using namespace seqan3::literals;

    auto char_to_char = seqan3::to_char('A');      // calls seqan3::custom::to_char('A')
    auto dna5_to_char = seqan3::to_char('A'_dna5); // calls .to_char() member
    std::cout << char_to_char << '\n';             // A
    std::cout << dna5_to_char << '\n';             // A
}
