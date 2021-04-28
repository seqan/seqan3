#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::cigar::operation letter{'M'_cigar_operation};

    letter.assign_char('D');
    seqan3::debug_stream << letter << '\n'; // prints "D"

    letter.assign_char('Z'); // Unknown characters are implicitly converted to M.
    seqan3::debug_stream << letter << '\n'; // prints "M"
}
