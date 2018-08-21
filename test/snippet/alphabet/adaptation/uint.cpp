#include <seqan3/alphabet/adaptation/uint.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

int main()
{

using namespace seqan3;
{
//! [assign_char]
auto l = assign_char(dna5{}, 'G');  // l is of type dna5
//! [assign_char]
(void) l;
}

{
//! [assign_rank]
auto l = assign_rank(dna5{}, 1);  // l is of type dna5 and == dna5::C
//! [assign_rank]
(void) l;
}

}
