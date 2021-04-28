#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::aa20 letter{'A'_aa20};

    letter.assign_char('C');
    seqan3::debug_stream << letter << '\n'; // prints "C"

    letter.assign_char('?'); // Unknown characters are implicitly converted to S.
    seqan3::debug_stream << letter << '\n'; // prints "S"
}
