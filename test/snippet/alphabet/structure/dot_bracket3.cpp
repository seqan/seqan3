#include <seqan3/alphabet/structure/dot_bracket3.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
//! [general]
// create vector
std::vector<dot_bracket3> vec{'.'_db3, ')'_db3, ')'_db3};
// modify and print
vec[1] = '('_db3;
for (dot_bracket3 chr : vec)
    debug_stream << to_char(chr);  // .()
debug_stream << "\n";
//! [general]

//! [string_literal]
std::vector<dot_bracket3> foo{".(..)."_db3};
std::vector<dot_bracket3> bar = ".(..)."_db3;
auto bax = ".(..)."_db3;
//! [string_literal]

//! [char_literal]
dot_bracket3 my_letter{'('_db3};

// does not work:
// dot_bracket3 my_letter{'('}; // <- char not implicitly convertible

my_letter.assign_char(')'); // <- assigns the char explicitly
//! [char_literal]
}
