#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/core/debug_stream.hpp>
int main()
{
    using seqan3::operator""_rna5;

    seqan3::rna5 my_letter{'A'_rna5};

    my_letter.assign_char('C');
    
    my_letter.assign_char('F'); // unknown characters are implicitly converted to N.
    if (my_letter.to_char() == 'N')
        seqan3::debug_stream << "yeah\n"; // "yeah";
}
