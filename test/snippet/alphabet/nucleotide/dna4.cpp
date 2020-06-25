#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_dna4;

    seqan3::dna4 my_letter{'C'_dna4};

    my_letter.assign_char('F'); // characters other than IUPAC characters are implicitly converted to A.
    seqan3::debug_stream << my_letter; // prints "F"

    // IUPAC characters are implicitly converted to their best fitting representative
    seqan3::debug_stream << my_letter.assign_char('R'); // prints "A"
    seqan3::debug_stream << my_letter.assign_char('Y'); // prints "C"
    seqan3::debug_stream << my_letter.assign_char('S'); // prints "C"
    seqan3::debug_stream << my_letter.assign_char('W'); // prints "A"
    seqan3::debug_stream << my_letter.assign_char('K'); // prints "G"
    seqan3::debug_stream << my_letter.assign_char('M'); // prints "A"
    seqan3::debug_stream << my_letter.assign_char('B'); // prints "C"
    seqan3::debug_stream << my_letter.assign_char('D'); // prints "A"
    seqan3::debug_stream << my_letter.assign_char('H'); // prints "A"
    seqan3::debug_stream << my_letter.assign_char('V'); // prints "A"

    my_letter.assign_char('a'); // lower case letters are the same as their upper case equivalent
    seqan3::debug_stream << my_letter; // prints "A";
}
