#include <seqan3/alphabet/aminoacid/aa10li.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/core/debug_stream.hpp>

using seqan3::operator""_aa10li;
using seqan3::operator""_aa27;

int main()
{
    // Construction of aa10li amino acids from character
    seqan3::aa10li my_letter{'A'_aa10li};

    my_letter.assign_char('C');
    my_letter.assign_char('?'); // all unknown characters are converted to 'A'_aa10li implicitly

    seqan3::debug_stream << my_letter.to_char() << '\n'; // A

    // Convert aa27 alphabet to aa10_murphy
    seqan3::aa27_vector v1{"ALRSTXOUMP"_aa27};
    auto v2 = v1 | std::views::transform([] (auto const & in) { return static_cast<seqan3::aa10li>(in); });

    seqan3::debug_stream << v2 << '\n'; // AJKAASKCJP

    return 0;
}
