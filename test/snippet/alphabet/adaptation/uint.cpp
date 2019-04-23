#include <seqan3/alphabet/adaptation/uint.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

using namespace seqan3;

int main()
{

{
//! [assign_char]
auto my_letter = assign_char_to('G', dna5{});  // my_letter is of type dna5 and == 'G'_dna5
//! [assign_char]
(void) my_letter;
}

{
//! [assign_rank]
auto my_letter = assign_rank_to(1, dna5{});  // my_letter is of type dna5 and == 'C'_dna5
//! [assign_rank]
(void) my_letter;
}

}
