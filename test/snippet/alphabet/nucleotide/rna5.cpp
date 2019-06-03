#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{

{
//! [code]
rna5 my_letter{'A'_rna5};
// doesn't work:
// rna5 my_letter{'A'};

my_letter.assign_char('C'); // <- this does!

my_letter.assign_char('F'); // unknown characters are implicitly converted to N.
if (my_letter.to_char() == 'N')
    debug_stream << "yeah\n"; // "yeah";
//! [code]
}

{
//! [operator""_rna5]
// these don't work:
// rna5_vector foo{"ACGTTA"};
// rna5_vector bar = "ACGTTA";

// but these do:
rna5_vector foo{"ACGTTA"_rna5};
rna5_vector bar = "ACGTTA"_rna5;
auto bax = "ACGTTA"_rna5;
//! [operator""_rna5]
}

}
