#include <seqan3/alphabet/structure/structured_aa.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>

int main()
{
using namespace seqan3;
{
//! [general]
structured_aa<aa27, dssp9> l{aa27::W, dssp9::B};
std::cout << int(to_rank(l)) << ' '
          << int(to_rank(get<0>(l))) << ' '
          << int(to_rank(get<1>(l))) << '\n';
// 49 22 1

std::cout << to_char(l) << ' '
          << to_char(get<0>(l)) << ' '
          << to_char(get<1>(l)) << '\n';
// W W B

// modify via structured bindings and references:
auto & [ seq_l, structure_l ] = l;
seq_l = aa27::V;
std::cout << to_char(l) << '\n';
// V
//! [general]
(void) structure_l;
}
}
