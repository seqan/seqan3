#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{

{
//! [construction]
aa27 my_letter{'A'_aa27};
// doesn't work:
// aa27 my_letter{'A'};

my_letter.assign_char('C'); // <- this does!

my_letter.assign_char('?'); // all unknown characters are converted to 'X'_aa27 implicitly
if (my_letter.to_char() == 'X')
    debug_stream << "yeah\n"; // "yeah";
//! [construction]
}

{
//! [char_literal]
// doesn't work:
// aa27 acid{'A'};

// this does:
aa27 acid1{'A'_aa27};
auto acid2 = 'Y'_aa27; // type = aa27
//! [char_literal]
}

{
//! [literal]
// these don't work:
// aa27_vector foo{"ABFUYR"};
// aa27_vector bar = "ABFUYR";

// but these do:
aa27_vector foo{"ABFUYR"_aa27};
aa27_vector bar = "ABFUYR"_aa27;
auto bax = "ABFUYR"_aa27;
//! [literal]
}

}
