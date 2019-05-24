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

using seqan3::operator""_wuss51;
//! [is_pair_open]
bool is_opening_char = '{'_wuss51.is_pair_open();       // true
is_opening_char = seqan3::is_pair_open('.'_wuss51);     // false
//! [is_pair_open]

//! [is_pair_close]
bool is_closing_char = '}'_wuss51.is_pair_close();      // true
is_closing_char = seqan3::is_pair_close('.'_wuss51);    // false
//! [is_pair_close]

//! [is_unpaired]
bool is_unpaired_char = '.'_wuss51.is_unpaired();       // true
is_unpaired_char = seqan3::is_unpaired('{'_wuss51);     // false
//! [is_unpaired]

//! [max_pseudoknot_depth]
uint8_t max_depth = wuss51::max_pseudoknot_depth;       // 22
max_depth = seqan3::max_pseudoknot_depth<wuss51>;       // 22
//! [max_pseudoknot_depth]

//! [pseudoknot_id]
auto pk_opt = '.'_wuss51.pseudoknot_id();               // std::optional -> false
pk_opt = seqan3::pseudoknot_id('{'_wuss51);             // std::optional -> true: 3

if (pk_opt)
    debug_stream << *pk_opt;                            // 3
//! [pseudoknot_id]

(void) is_opening_char;
(void) is_closing_char;
(void) is_unpaired_char;
(void) max_depth;
}
