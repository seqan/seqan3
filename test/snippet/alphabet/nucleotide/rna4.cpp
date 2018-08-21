#include <seqan3/alphabet/nucleotide/rna4.hpp>

int main()
{
    using namespace seqan3;

//! [code]
rna4 my_letter{rna4::A};
// doesn't work:
// rna4 my_letter{'A'};

my_letter.assign_char('C'); // <- this does!

my_letter.assign_char('F'); // converted to A internally
if (my_letter.to_char() == 'A')
    std::cout << "yeah\n"; // "yeah";
//! [code]

//! [operator""_rna4]
// these don't work:
// rna4_vector foo{"ACGTTA"};
// rna4_vector bar = "ACGTTA";

// but these do:
using namespace seqan3::literal;
rna4_vector foo{"ACGTTA"_rna4};
rna4_vector bar = "ACGTTA"_rna4;
auto bax = "ACGTTA"_rna4;
//! [operator""_rna4]
}
