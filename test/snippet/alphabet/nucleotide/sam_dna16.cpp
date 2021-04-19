#include <seqan3/alphabet/nucleotide/sam_dna16.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_dna16sam;

    seqan3::dna16sam my_letter{'A'_dna16sam};

    my_letter.assign_char('=');

    my_letter.assign_char('F'); // unknown characters are implicitly converted to N.
    seqan3::debug_stream << my_letter << '\n'; // "N";
}
