#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::cigar letter{12, 'M'_cigar_operation};

    letter.assign_string("10D");
    seqan3::debug_stream << letter << '\n'; // prints "10D"

    letter.assign_string("20Z"); // Unknown strings are implicitly converted to 0P.
    seqan3::debug_stream << letter << '\n'; // prints "0P"
}
