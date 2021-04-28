#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna15 letter{'A'_dna15};

    letter.assign_char('C');
    seqan3::debug_stream << letter << '\n'; // prints "C"

    letter.assign_char('F'); // Unknown characters are implicitly converted to N.
    seqan3::debug_stream << letter << '\n'; // prints "N"
}
