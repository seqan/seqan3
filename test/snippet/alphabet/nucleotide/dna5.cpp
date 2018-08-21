#include <seqan3/alphabet/nucleotide/dna5.hpp>

using namespace seqan3;

int main()
{

{
//! [code]
dna5 my_letter{dna5::A};
// doesn't work:
// dna5 my_letter{'A'};

my_letter.assign_char('C'); // <- this does!

my_letter.assign_char('F'); // converted to N internally
if (my_letter.to_char() == 'N')
    std::cout << "yeah\n"; // "yeah";
//! [code]
}

{
//! [operator""_dna5]
// these don't work:
// dna5_vector foo{"ACGTTA"};
// dna5_vector bar = "ACGTTA";

// but these do:
using namespace seqan3::literal;
dna5_vector foo{"ACGTTA"_dna5};
dna5_vector bar = "ACGTTA"_dna5;
auto bax = "ACGTTA"_dna5;
//! [operator""_dna5]
}

}
