#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{

{
//! [ambiguity]
// does not work:
// dna4 my_letter{0};      // we want to set the default, an A
// dna4 my_letter{'A'};    // we also want to set an A, but we are setting value 65

// std::cout << my_letter; // you expect 'A', but how would you access the number?
//! [ambiguity]
}

{
//! [nonambiguous]
dna4 my_letter;
assign_rank_to(0, my_letter);       // assign an A via rank interface
assign_char_to('A', my_letter);     // assign an A via char interface
//! [nonambiguous]

//! [print]
std::cout << to_char(my_letter);            // prints 'A'
std::cout << (unsigned)to_rank(my_letter);  // prints 0
// we have to add the cast here, because uint8_t is also treated as a char type by default :(

// Using SeqAn's debug_stream:
debug_stream << to_char(my_letter);         // prints 'A'
debug_stream << my_letter;                  // prints 'A' (calls to_char() automatically!)
debug_stream << to_rank(my_letter);         // prints 0   (casts uint8_t to unsigned automatically!)
//! [print]
}

{
//! [literal]
dna4        my_letter = 'A'_dna4;           // identical to assign_char_to('A', my_letter);
dna4_vector my_seq    = "ACGT"_dna4;        // identical to calling assign_char for each element
//! [literal]
}

}
