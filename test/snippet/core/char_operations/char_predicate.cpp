#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/char_operations/predicate.hpp>

using namespace seqan3;

int main()
{
{
//! [is_in_interval]
is_in_interval<'A', 'G'>('C'); // returns true

auto constexpr my_check = is_in_interval<'A', 'G'>;
my_check('H');  // returns false
//! [is_in_interval]
}

{
//! [is_in_alphabet]
is_in_alphabet<dna4>('C');  // returns true

auto constexpr my_check = is_in_alphabet<dna4>;
my_check('U');  // returns false, because it comes out as 'T'
//! [is_in_alphabet]
}

{
//! [is_char]
is_char<'C'>('C');  // returns true

auto constexpr my_check = is_char<'C'>;
my_check('c');  // returns false, because case is different
//! [is_char]
}

//! [is_eof]
is_eof(EOF);  // returns true
is_eof('C');  // returns false
//! [is_eof]

//! [is_cntrl]
is_cntrl('\0');  // returns true.
//! [is_cntrl]

//! [is_print]
is_print(' ');  // returns true.
//! [is_print]

//! [is_space]
is_space('\n');  // returns true.
//! [is_space]

//! [is_blank]
is_blank('\t');  // returns true.
//! [is_blank]

//! [is_graph]
is_graph('%');  // returns true.
//! [is_graph]

//! [is_punct]
is_punct(':');  // returns true.
//! [is_punct]

//! [is_alnum]
is_alnum('9');  // returns true.
//! [is_alnum]

//! [is_alpha]
is_alpha('z');  // returns true.
//! [is_alpha]

//! [is_upper]
is_upper('K');  // returns true.
//! [is_upper]

//! [is_lower]
is_lower('a');  // returns true.
//! [is_lower]

//! [is_digit]
is_digit('1');  // returns true.
//! [is_digit]

//! [is_xdigit]
is_xdigit('e');  // returns true.
//! [is_xdigit]

{
//! [operator]
char chr{'1'};
auto constexpr my_cond = is_char<'%'> || is_digit;
bool is_percent = my_cond(chr); // is_percent == true
//! [operator]
(void) is_percent;
}

}
