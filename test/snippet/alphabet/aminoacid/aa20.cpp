#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{

{
//! [construction]
aa20 my_letter{'A'_aa20};
// doesn't work:
// aa20 my_letter{'A'};

my_letter.assign_char('C'); // <- this does!

my_letter.assign_char('?'); // all unknown characters are converted to 'A'_aa20 implicitly

if (my_letter.to_char() == 'A')
    debug_stream << "yeah\n"; // "yeah";
//! [construction]
}

{
//! [char_literal]
// doesn't work:
// aa20 acid{'A'};

// this does:
aa20 acid1{'A'_aa20};
auto acid2 = 'Y'_aa20; // type = aa20
//! [char_literal]
}

{
//! [literal]
// these don't work:
// aa20_vector foo{"ABFUYR"};
// aa20_vector bar = "ABFUYR";

// but these do:
aa20_vector foo{"ABFUYR"_aa20};
aa20_vector bar = "ABFUYR"_aa20;
auto bax = "ABFUYR"_aa20;
//! [literal]
}

}
