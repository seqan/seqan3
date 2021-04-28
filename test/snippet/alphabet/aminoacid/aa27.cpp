#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::aa27 letter{'A'_aa27};

    letter.assign_char('C');
    seqan3::debug_stream << letter << '\n'; // prints "C"

    letter.assign_char('?'); // Unknown characters are implicitly converted to X.
    seqan3::debug_stream << letter << '\n'; // prints "X"
}
