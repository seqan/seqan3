#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_rna15;

    seqan3::rna15 my_letter{'A'_rna15};

    my_letter.assign_char('C');

    my_letter.assign_char('F'); // unknown characters are implicitly converted to N.
    if (my_letter.to_char() == 'N')
        seqan3::debug_stream << "yeah\n"; // "yeah";
}
