#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::wuss51 letter{':'_wuss51};

    letter.assign_char('~');
    seqan3::debug_stream << letter << '\n'; // prints "~"

    letter.assign_char('#'); // Unknown characters are implicitly converted to ';'.
    seqan3::debug_stream << letter << '\n'; // prints ";"
}
