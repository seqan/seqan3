#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

using namespace seqan3;

int main()
{

{
//! [ambiguity]
// does not work:
// dna4 my_letter{0};      // we want to set the default, an A
// dna4 my_letter{'A'};    // we also want to set an A, but we are setting value 65

// debug_stream << my_letter; // you expect 'A', but how would you access the number?
//! [ambiguity]
}

{
//! [nonambiguous]
dna4 my_letter;
assign_rank(my_letter, 0);       // assign an A via rank interface
assign_char(my_letter, 'A');     // assign an A via char interface
my_letter = 'A'_dna4;             // some alphabets (BUT NOT ALL!) also provide an enum-like interface

debug_stream << to_char(my_letter);            // prints 'A'
debug_stream << (unsigned)to_rank(my_letter);  // prints 0
// we have to add the cast here, because uint8_t is also treated as a char type by default :(
//! [nonambiguous]
}

}
