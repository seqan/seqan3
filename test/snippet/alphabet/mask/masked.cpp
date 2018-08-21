#include <seqan3/alphabet/mask/masked.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <iostream>

int main()
{
using namespace seqan3;
{
//! [general]
masked<dna4> dna4_masked{};
masked<dna4> dna4_another_masked{dna4::A, mask::UNMASKED}; // create a dna4 masked alphabet with an unmasked A

dna4_masked.assign_char('a'); // assigns a masked dna4::A

if (dna4_masked.to_char() != dna4_another_masked.to_char())
    std::cout << dna4_masked.to_char() << " is not the same as " << dna4_another_masked.to_char() << "\n";
//! [general]
}
}
