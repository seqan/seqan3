#include <seqan3/alphabet/aminoacid/aa27.hpp>

using namespace seqan3;

int main()
{

{
//! [construction]
aa27 my_letter{aa27::A};
// doesn't work:
// aa27 my_letter{'A'};

my_letter.assign_char('C'); // <- this does!

my_letter.assign_char('?'); // converted to X internally
if (my_letter.to_char() == 'X')
    std::cout << "yeah\n"; // "yeah";
//! [construction]
}

{
//! [literal]
// these don't work:
// aa27_vector foo{"ABFUYR"};
// aa27_vector bar = "ABFUYR";

// but these do:
using namespace seqan3::literal;
aa27_vector foo{"ABFUYR"_aa27};
aa27_vector bar = "ABFUYR"_aa27;
auto bax = "ABFUYR"_aa27;
//! [literal]
}

}
