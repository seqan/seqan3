#include <seqan3/alphabet/structure/dssp9.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dssp9 letter{'H'_dssp9};

    letter.assign_char('B');
    seqan3::debug_stream << letter << '\n'; // prints "B"

    letter.assign_char('F'); // Unknown characters are implicitly converted to 'X'.
    seqan3::debug_stream << letter << '\n'; // prints "X"
}
