#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

using namespace seqan3;

int main()
{

{
//! [code]
dna15 my_letter{dna15::A};
// doesn't work:
// dna15 my_letter{'A'};

my_letter.assign_char('C'); // <- this does!

my_letter.assign_char('F'); // unknown characters are implicitly converted to N.
if (my_letter.to_char() == 'N')
    debug_stream << "yeah\n"; // "yeah";
//! [code]
}

{
//! [operator""_dna15]
// these don't work:
// dna15_vector foo{"ACGTTA"};
// dna15_vector bar = "ACGTTA";

// but these do:
using namespace seqan3::literal;
dna15_vector foo{"ACGTTA"_dna15};
dna15_vector bar = "ACGTTA"_dna15;
auto bax = "ACGTTA"_dna15;
//! [operator""_dna15]
}

}
