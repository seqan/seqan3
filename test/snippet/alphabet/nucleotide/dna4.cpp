#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_dna4;

    seqan3::dna4 my_letter{'A'_dna4};

    my_letter.assign_char('C');

    my_letter.assign_char('F'); // unknown characters are implicitly converted to A.
    if (my_letter.to_char() == 'A')
        seqan3::debug_stream << "yeah\n"; // "yeah";

}
