#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_aa27;

    seqan3::aa27 my_letter{'A'_aa27};

    my_letter.assign_char('C');
    
    my_letter.assign_char('?'); // all unknown characters are converted to 'X'_aa27 implicitly
    if (my_letter.to_char() == 'X')
        seqan3::debug_stream << "yeah\n"; // "yeah";
}
