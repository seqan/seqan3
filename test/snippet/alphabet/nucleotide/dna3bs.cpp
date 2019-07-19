#include <seqan3/alphabet/nucleotide/dna3bs.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_dna3bs;

    seqan3::dna3bs my_letter{'A'_dna3bs};

    my_letter.assign_char('C'); // all C will be converted to T.
    if (my_letter.to_char() == 'T')
        seqan3::debug_stream << "yeah\n"; // "yeah";

    my_letter.assign_char('F'); // unknown characters are implicitly converted to A.
    if (my_letter.to_char() == 'A')
        seqan3::debug_stream << "yeah\n"; // "yeah";

    return 0;
}
