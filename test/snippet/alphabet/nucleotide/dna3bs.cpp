#include <seqan3/alphabet/nucleotide/dna3bs.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna3bs letter{'A'_dna3bs};

    letter.assign_char('C'); // All C will be converted to T.
    seqan3::debug_stream << letter << '\n'; // prints "T"

    letter.assign_char('F'); // Unknown characters are implicitly converted to A.
    seqan3::debug_stream << letter << '\n'; // prints "A"
}
