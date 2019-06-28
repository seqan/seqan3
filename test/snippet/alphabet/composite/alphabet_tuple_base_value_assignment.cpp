#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>

int main()
{
    using seqan3::operator""_dna4;

    seqan3::qualified<seqan3::dna4, seqan3::phred42> letter1{'T'_dna4, seqan3::phred42{1}};

    letter1 = 'C'_dna4;           // yields {'C'_dna4, seqan3::phred42{1}}
    letter1 = seqan3::phred42{2}; // yields {'C'_dna4, seqan3::phred42{2}}
}
