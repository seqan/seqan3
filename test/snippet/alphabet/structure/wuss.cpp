#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
//! [general]
// create vector
std::vector<wuss51> vec{'.'_wuss51, '>'_wuss51, '>'_wuss51};
// modify and print
vec[1] = '<'_wuss51;
for (wuss51 chr : vec)
    debug_stream << to_char(chr);  // .<>
debug_stream << "\n";
//! [general]

//! [string_literal]
std::vector<wuss<>> foo{".<..>."_wuss51};
std::vector<wuss<>> bar = ".<..>."_wuss51;
auto bax = ".<..>."_wuss51;
//! [string_literal]

//! [char_literal]
wuss51 my_letter{'~'_wuss51};

// does not work:
// wuss51 my_letter{'~'}; // <- char not implicitly convertible

// works for each wuss alphabet size:
my_letter.assign_char('<'); // <- assigns the char explicitly
//! [char_literal]
}
