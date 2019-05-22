#include <seqan3/alphabet/structure/structured_rna.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/structure/dot_bracket3.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
//! [general]
structured_rna<rna4, dot_bracket3> letter{'G'_rna4, '('_db3};
debug_stream << int(to_rank(letter)) << ' '
          << int(to_rank(get<0>(letter))) << ' '
          << int(to_rank(get<1>(letter))) << '\n';
// 6 2 1

debug_stream << to_char(letter) << ' '
          << to_char(get<0>(letter)) << ' '
          << to_char(get<1>(letter)) << '\n';
// G G (

// modify via structured bindings and references:
auto & [ seq_l, structure_l ] = letter;
seq_l = 'U'_rna4;
debug_stream << to_char(letter) << '\n';
// U
//! [general]
(void) structure_l;
}
