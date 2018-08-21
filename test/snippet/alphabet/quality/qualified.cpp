#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3;

int main()
{
//! [general]
qualified<dna4, phred42> l{dna4::A, 7};
std::cout << int(to_rank(l)) << ' '
          << int(to_rank(get<0>(l))) << ' '
          << int(to_rank(get<1>(l))) << '\n';
// 28 0 7

std::cout << to_char(l) << ' '
          << to_char(get<0>(l)) << ' '
          << to_char(get<1>(l)) << '\n';
// A A (

std::cout << int(to_phred(l)) << ' '
//        << int(to_phred(get<0>(l))) << ' ' // dna4 doesn't have a phred
          << int(to_phred(get<1>(l))) << '\n';
// 7 7

// modify via structured bindings and references:
auto & [ seq_l, qual_l ] = l;
seq_l = dna4::G;
std::cout << to_char(l) << '\n';
// G
//! [general]
(void) qual_l;
}
