#include <seqan3/alphabet/aminoacid/aa10murphy.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/convert.hpp>

using seqan3::operator""_aa10murphy;
using seqan3::operator""_aa27;

int main()
{
    // Construction of aa10murphy amino acids from character
    seqan3::aa10murphy my_letter{'A'_aa10murphy};

    my_letter.assign_char('C');
    my_letter.assign_char('?'); // all unknown characters are converted to 'S'_aa10murphy implicitly

    if (my_letter.to_char() == 'S')
        seqan3::debug_stream << "yeah\n"; // "yeah";

    // Convert aa27 alphabet to aa10_murphy
    seqan3::aa27_vector v1{"AVRSTXOUB"_aa27};
    auto v2 = v1 | seqan3::view::convert<seqan3::aa10murphy>; // AIKSSSKCB

    seqan3::debug_stream << v2 << "\n";

    return 0;
}
