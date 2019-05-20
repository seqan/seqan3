#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{

{
//! [code]
rna4 my_letter{'A'_rna4};
// doesn't work:
// rna4 my_letter{'A'};

my_letter.assign_char('C'); // <- this does!

my_letter.assign_char('F'); // unknown characters are implicitly converted to A.
if (my_letter.to_char() == 'A')
    debug_stream << "yeah\n"; // "yeah";
//! [code]
}

{
//! [operator""_rna4]
// these don't work:
// rna4_vector foo{"ACGTTA"};
// rna4_vector bar = "ACGTTA";

// but these do:
rna4_vector foo{"ACGTTA"_rna4};
rna4_vector bar = "ACGTTA"_rna4;
auto bax = "ACGTTA"_rna4;
//! [operator""_rna4]
}

}
