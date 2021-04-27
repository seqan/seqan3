#include <seqan3/alphabet/nucleotide/dna16sam.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna16sam letter{'A'_dna16sam};

    letter.assign_char('=');
    seqan3::debug_stream << letter << '\n'; // prints "="

    letter.assign_char('F'); // Unknown characters are implicitly converted to N.
    seqan3::debug_stream << letter << '\n'; // "N";
}
