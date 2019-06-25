#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_dna4;

    seqan3::qualified<seqan3::dna4, seqan3::phred42> letter1{'C'_dna4};
    // creates {'C'_dna4, seqan3::phred42{0}}
    seqan3::qualified<seqan3::dna4, seqan3::phred42> letter2{seqan3::phred42{1}};
    // creates {'A'_dna4, seqan3::phred42{1}}

    if (letter1 == letter2)
        seqan3::debug_stream << "yeah\n"; // yeah
}
