#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
//! [general]
qualified<dna4, phred42> letter{'A'_dna4, phred42{7}};
debug_stream << int(to_rank(letter)) << ' '
          << int(to_rank(get<0>(letter))) << ' '
          << int(to_rank(get<1>(letter))) << '\n';
// 28 0 7

debug_stream << to_char(letter) << ' '
          << to_char(get<0>(letter)) << ' '
          << to_char(get<1>(letter)) << '\n';
// A A (

debug_stream << int(to_phred(letter)) << ' '
//        << int(to_phred(get<0>(letter))) << ' ' // dna4 doesn't have a phred
          << int(to_phred(get<1>(letter))) << '\n';
// 7 7

// modify via structured bindings and references:
auto & [ seq_l, qual_l ] = letter;
seq_l = 'G'_dna4;
debug_stream << to_char(letter) << '\n';
// G
//! [general]
(void) qual_l;
}
