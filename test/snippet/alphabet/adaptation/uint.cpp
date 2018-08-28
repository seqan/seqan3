#include <seqan3/alphabet/adaptation/uint.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

using namespace seqan3;

int main()
{

{
//! [assign_char]
auto my_letter = assign_char(dna5{}, 'G');  // my_letter is of type dna5 and == dna5::G
//! [assign_char]
(void) my_letter;
}

{
//! [assign_rank]
auto my_letter = assign_rank(dna5{}, 1);  // my_letter is of type dna5 and == dna5::C
//! [assign_rank]
(void) my_letter;
}

}
