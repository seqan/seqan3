#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

int main()
{
    using namespace seqan3::literals;

    auto r2 = seqan3::to_char('A');         // calls seqan3::custom::to_char('A'); r2 == 'A'
    auto r3 = seqan3::to_char('A'_dna5);    // calls .to_char() member; r3 == 'A'
}
