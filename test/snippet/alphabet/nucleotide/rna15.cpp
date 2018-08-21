#include <seqan3/alphabet/nucleotide/rna15.hpp>

using namespace seqan3;

int main()
{

{
//! [code]
rna15 my_letter{rna15::A};
// doesn't work:
// rna15 my_letter{'A'};

my_letter.assign_char('C'); // <- this does!

my_letter.assign_char('F'); // converted to A internally
if (my_letter.to_char() == 'N')
    std::cout << "yeah\n"; // "yeah";
//! [code]
}

{
//! [operator""_rna15]
// these don't work:
// rna15_vector foo{"ACGTTA"};
// rna15_vector bar = "ACGTTA";

// but these do:
using namespace seqan3::literal;
rna15_vector foo{"ACGTTA"_rna15};
rna15_vector bar = "ACGTTA"_rna15;
auto bax = "ACGTTA"_rna15;
//! [operator""_rna15]
}

}
