#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/concept.hpp>

int main()
{
    char c = '!';
    seqan3::assign_char_to('?', c);     // calls seqan3::adaptation::assign_char_to('A', c)

    seqan3::dna5 d{};
    seqan3::assign_char_to('A', d);     // calls .assign_char('A') member

    // also works for temporaries:
    seqan3::dna5 d2 = seqan3::assign_char_to('A', seqan3::dna5{});

    // invalid/unknown characters are converted:
    seqan3::dna5 d3 = seqan3::assign_char_to('!', seqan3::dna5{}); // == 'N'_dna5
}
