#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

int main()
{
    using seqan3::operator""_dna4;

    seqan3::gapped<seqan3::dna4> gapped_letter{};
    seqan3::gapped<seqan3::dna4> converted_letter{'C'_dna4};
    seqan3::gapped<seqan3::dna4> gap_letter{seqan3::gap{}};

    seqan3::gapped<seqan3::dna4>{}.assign_char('C');
    seqan3::gapped<seqan3::dna4>{}.assign_char('-'); // gap character
    seqan3::gapped<seqan3::dna4>{}.assign_char('K'); // unknown characters map to the default/unknown
                                                     // character of the given alphabet type (i.e. A of dna4)
}
