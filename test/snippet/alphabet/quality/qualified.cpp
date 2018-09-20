#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3;

int main()
{
//! [general]
qualified<dna4, phred42> letter{dna4::A, phred42{7}};
std::cout << int(to_rank(letter)) << ' '
          << int(to_rank(get<0>(letter))) << ' '
          << int(to_rank(get<1>(letter))) << '\n';
// 28 0 7

std::cout << to_char(letter) << ' '
          << to_char(get<0>(letter)) << ' '
          << to_char(get<1>(letter)) << '\n';
// A A (

std::cout << int(to_phred(letter)) << ' '
//        << int(to_phred(get<0>(letter))) << ' ' // dna4 doesn't have a phred
          << int(to_phred(get<1>(letter))) << '\n';
// 7 7

// modify via structured bindings and references:
auto & [ seq_l, qual_l ] = letter;
seq_l = dna4::G;
std::cout << to_char(letter) << '\n';
// G
//! [general]
(void) qual_l;
}
