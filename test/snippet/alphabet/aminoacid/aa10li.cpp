#include <seqan3/alphabet/aminoacid/aa10li.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/convert.hpp>

using seqan3::operator""_aa10li;
using seqan3::operator""_aa27;

int main()
{
    // Construction of aa10li amino acids from character
    seqan3::aa10li my_letter{'A'_aa10li};

    my_letter.assign_char('C');
    my_letter.assign_char('?'); // all unknown characters are converted to 'S'_aa10li implicitly

    if (my_letter.to_char() == 'S')
        seqan3::debug_stream << "yeah\n"; // "yeah";

    // Convert aa27 alphabet to aa10_murphy
    seqan3::aa27_vector v1{"ALRSTXOUMP"_aa27};
    auto v2 = v1 | seqan3::view::convert<seqan3::aa10li>; // AJKAASKCJP

    seqan3::debug_stream << v2 << "\n";

    return 0;
}
