#include <seqan3/alphabet/structure/structured_rna.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/structure/dot_bracket3.hpp>

int main()
{
using namespace seqan3;
{
//! [general]
structured_rna<rna4, dot_bracket3> l{rna4::G, dot_bracket3::PAIR_OPEN};
std::cout << int(to_rank(l)) << ' '
          << int(to_rank(get<0>(l))) << ' '
          << int(to_rank(get<1>(l))) << '\n';
// 6 2 1

std::cout << to_char(l) << ' '
          << to_char(get<0>(l)) << ' '
          << to_char(get<1>(l)) << '\n';
// G G (

// modify via structured bindings and references:
auto & [ seq_l, structure_l ] = l;
seq_l = rna4::U;
std::cout << to_char(l) << '\n';
// U
//! [general]
(void) structure_l;
}
}
