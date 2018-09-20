#include <seqan3/alphabet/structure/structured_aa.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>

using namespace seqan3;
int main()
{
//! [general]
structured_aa<aa27, dssp9> letter{aa27::W, dssp9::B};
std::cout << int(to_rank(letter)) << ' '
          << int(to_rank(get<0>(letter))) << ' '
          << int(to_rank(get<1>(letter))) << '\n';
// 49 22 1

std::cout << to_char(letter) << ' '
          << to_char(get<0>(letter)) << ' '
          << to_char(get<1>(letter)) << '\n';
// W W B

// modify via structured bindings and references:
auto & [ seq_l, structure_l ] = letter;
seq_l = aa27::V;
std::cout << to_char(letter) << '\n';
// V
//! [general]
(void) structure_l;
}
