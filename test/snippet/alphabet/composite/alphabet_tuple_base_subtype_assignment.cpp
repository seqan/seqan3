#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>

int main()
{
    using seqan3::operator""_dna4;
    using seqan3::operator""_rna4;

    seqan3::qualified<seqan3::dna4, seqan3::phred42> letter1{'T'_dna4, seqan3::phred42{1}};

    letter1 = 'C'_rna4; // yields {'C'_dna4, seqan3::phred42{1}}
}
