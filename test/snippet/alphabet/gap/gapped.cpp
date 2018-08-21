#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3;

int main()
{
//! [general]
gapped<dna4> gapped_letter{};
gapped<dna4> converted_letter{dna4::C};
// doesn't work:
// gapped<dna4> my_letter{'A'};

gapped<dna4>{}.assign_char('C'); // <- this does!
gapped<dna4>{}.assign_char('-'); // gap character
gapped<dna4>{}.assign_char('K'); // unknown characters map to the default/unknown
                                 // character of the given alphabet type (i.e. A of dna4)
//! [general]
}
