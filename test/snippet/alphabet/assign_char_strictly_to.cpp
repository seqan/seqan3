#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

int main()
{
    char c = '!';
    seqan3::assign_char_strictly_to('?', c);     // calls seqan3::custom::assign_char_strictly_to('A', c)

    seqan3::dna5 d{};
    seqan3::assign_char_strictly_to('A', d);     // calls .assign_char('A') member

    // also works for temporaries:
    seqan3::dna5 d2 = seqan3::assign_char_strictly_to('A', seqan3::dna5{});
}
