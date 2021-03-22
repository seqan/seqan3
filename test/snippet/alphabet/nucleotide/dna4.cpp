#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_dna4;

    seqan3::dna4 my_letter{'C'_dna4};

    my_letter.assign_char('F'); // characters other than IUPAC characters are implicitly converted to A.
    seqan3::debug_stream << my_letter << '\n'; // prints "A"

    // IUPAC characters are implicitly converted to their best fitting representative
    seqan3::debug_stream << my_letter.assign_char('R') << '\n'; // prints "A"
    seqan3::debug_stream << my_letter.assign_char('Y') << '\n'; // prints "C"
    seqan3::debug_stream << my_letter.assign_char('S') << '\n'; // prints "C"
    seqan3::debug_stream << my_letter.assign_char('W') << '\n'; // prints "A"
    seqan3::debug_stream << my_letter.assign_char('K') << '\n'; // prints "G"
    seqan3::debug_stream << my_letter.assign_char('M') << '\n'; // prints "A"
    seqan3::debug_stream << my_letter.assign_char('B') << '\n'; // prints "C"
    seqan3::debug_stream << my_letter.assign_char('D') << '\n'; // prints "A"
    seqan3::debug_stream << my_letter.assign_char('H') << '\n'; // prints "A"
    seqan3::debug_stream << my_letter.assign_char('V') << '\n'; // prints "A"

    my_letter.assign_char('a'); // lower case letters are the same as their upper case equivalent
    seqan3::debug_stream << my_letter << '\n'; // prints "A";
}
