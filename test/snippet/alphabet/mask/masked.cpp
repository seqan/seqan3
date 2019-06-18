#include <iostream>
#include <seqan3/alphabet/mask/masked.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
//! [general]
masked<dna4> dna4_masked{};
masked<dna4> dna4_another_masked{'A'_dna4, mask::UNMASKED}; // create a dna4 masked alphabet with an unmasked A

dna4_masked.assign_char('a'); // assigns a masked 'A'_dna4

if (dna4_masked.to_char() != dna4_another_masked.to_char())
    debug_stream << dna4_masked.to_char() << " is not the same as " << dna4_another_masked.to_char() << "\n";
//! [general]
}
