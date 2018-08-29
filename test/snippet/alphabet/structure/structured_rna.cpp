#include <seqan3/alphabet/structure/structured_rna.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/structure/dot_bracket3.hpp>

using namespace seqan3;

int main()
{
//! [general]
structured_rna<rna4, dot_bracket3> letter{rna4::G, dot_bracket3::PAIR_OPEN};
std::cout << int(to_rank(letter)) << ' '
          << int(to_rank(get<0>(letter))) << ' '
          << int(to_rank(get<1>(letter))) << '\n';
// 6 2 1

std::cout << to_char(letter) << ' '
          << to_char(get<0>(letter)) << ' '
          << to_char(get<1>(letter)) << '\n';
// G G (

// modify via structured bindings and references:
auto & [ seq_l, structure_l ] = letter;
seq_l = rna4::U;
std::cout << to_char(letter) << '\n';
// U
//! [general]
(void) structure_l;
}
