#include <seqan3/alphabet/nucleotide/dna4.hpp>

int main()
{
    using namespace seqan3;

//! [code]
dna4 my_letter{dna4::A};
// doesn't work:
// dna4 my_letter{'A'};

my_letter.assign_char('C'); // <- this does!

my_letter.assign_char('F'); // converted to A internally
if (my_letter.to_char() == 'A')
    std::cout << "yeah\n"; // "yeah";
//! [code]

//! [operator""_dna4]
// these don't work:
// dna4_vector foo{"ACGTTA"};
// dna4_vector bar = "ACGTTA";

// but these do:
using namespace seqan3::literal;
dna4_vector foo{"ACGTTA"_dna4};
dna4_vector bar = "ACGTTA"_dna4;
auto bax = "ACGTTA"_dna4;
//! [operator""_dna4]
}
