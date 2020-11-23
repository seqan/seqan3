#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>

int main()
{
//! [is_eof]
seqan3::is_eof(EOF);  // returns true
seqan3::is_eof('C');  // returns false
//! [is_eof]

//! [is_cntrl]
seqan3::is_cntrl('\0');  // returns true.
//! [is_cntrl]

//! [is_print]
seqan3::is_print(' ');  // returns true.
//! [is_print]

//! [is_space]
seqan3::is_space('\n');  // returns true.
//! [is_space]

//! [is_blank]
seqan3::is_blank('\t');  // returns true.
//! [is_blank]

//! [is_graph]
seqan3::is_graph('%');  // returns true.
//! [is_graph]

//! [is_punct]
seqan3::is_punct(':');  // returns true.
//! [is_punct]

//! [is_alnum]
seqan3::is_alnum('9');  // returns true.
//! [is_alnum]

//! [is_alpha]
seqan3::is_alpha('z');  // returns true.
//! [is_alpha]

//! [is_upper]
seqan3::is_upper('K');  // returns true.
//! [is_upper]

//! [is_lower]
seqan3::is_lower('a');  // returns true.
//! [is_lower]

//! [is_digit]
seqan3::is_digit('1');  // returns true.
//! [is_digit]

//! [is_xdigit]
seqan3::is_xdigit('e');  // returns true.
//! [is_xdigit]
}
