#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{

{
//! [code]
rna15 my_letter{'A'_rna15};
// doesn't work:
// rna15 my_letter{'A'};

my_letter.assign_char('C'); // <- this does!

my_letter.assign_char('F'); // unknown characters are implicitly converted to N.
if (my_letter.to_char() == 'N')
    debug_stream << "yeah\n"; // "yeah";
//! [code]
}

{
//! [operator""_rna15]
// these don't work:
// rna15_vector foo{"ACGTTA"};
// rna15_vector bar = "ACGTTA";

// but these do:
rna15_vector foo{"ACGTTA"_rna15};
rna15_vector bar = "ACGTTA"_rna15;
auto bax = "ACGTTA"_rna15;
//! [operator""_rna15]
}

}
