#include <seqan3/alphabet/structure/dssp9.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
//! [general]
// create vector
std::vector<dssp9> vec{'E'_dssp9, 'H'_dssp9, 'H'_dssp9, 'H'_dssp9, 'T'_dssp9, 'G'_dssp9};
// modify and print
vec[1] = 'C'_dssp9;
for (dssp9 chr : vec)
    debug_stream << to_char(chr);  // ECHHTG
debug_stream << "\n";
//! [general]

//! [string_literal]
std::vector<dssp9> foo{"EHHHHT"_dssp9};
std::vector<dssp9> bar = "EHHHHT"_dssp9;
auto bax = "EHHHHT"_dssp9;
//! [string_literal]


//! [char_literal]
dssp9 my_letter{'I'_dssp9};

// does not work:
// dssp9 my_letter{'I'}; // <- char not implicitly convertible

my_letter.assign_char('G'); // <- assigns the char explicitly
//! [char_literal]
}
