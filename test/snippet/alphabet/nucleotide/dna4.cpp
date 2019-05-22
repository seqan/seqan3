#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{

{
//! [code]
dna4 my_letter{'A'_dna4};
// doesn't work:
// dna4 my_letter{'A'};

my_letter.assign_char('C'); // <- this does!

my_letter.assign_char('F'); // unknown characters are implicitly converted to A.
if (my_letter.to_char() == 'A')
    debug_stream << "yeah\n"; // "yeah";
//! [code]
}

{
//! [operator""_dna4]
// these don't work:
// dna4_vector foo{"ACGTTA"};
// dna4_vector bar = "ACGTTA";

// but these do:
dna4_vector foo{"ACGTTA"_dna4};
dna4_vector bar = "ACGTTA"_dna4;
auto bax = "ACGTTA"_dna4;
//! [operator""_dna4]
}

}
