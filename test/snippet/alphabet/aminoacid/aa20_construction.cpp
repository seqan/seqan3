#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_aa20;

    seqan3::aa20 my_letter{'A'_aa20};
    
    my_letter.assign_char('C');

    my_letter.assign_char('?'); // all unknown characters are converted to 'A'_aa20 implicitly

    if (my_letter.to_char() == 'A')
        seqan3::debug_stream << "yeah\n"; // "yeah";
}
