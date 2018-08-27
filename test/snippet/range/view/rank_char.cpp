#include <seqan3/range/view/char_to.hpp>
#include <seqan3/range/view/rank_to.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/range/view/to_rank.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>

#include <vector>


using namespace seqan3;
using namespace seqan3::literal;

int main()
{

{
//! [char_to]
std::string s{"ACTTTGATAN"};
auto v1 = s | view::char_to<dna4>; // == "ACTTTGATAA"_dna4
auto v2 = s | view::char_to<dna5>; // == "ACTTTGATAN"_dna5
//! [char_to]
(void) v1;
(void) v2;
}

{
//! [rank_to]
std::vector<int> vec{0, 1, 3, 3, 3, 2, 0, 3, 0};
auto v1 = vec | view::rank_to<dna4>; // == "ACTTTGATA"_dna4
auto v2 = vec | view::rank_to<dna5>; // == "ACTTTGATA"_dna5
//! [rank_to]
(void) v1;
(void) v2;
}

{
//! [to_char]
dna4_vector vec = "ACTTTGATA"_dna4;
auto v = vec | view::to_char;
std::cout << v << '\n'; // [A,C,T,T,T,G,A,T,A]

std::vector<phred42> qvec{{0}, {7}, {5}, {3}, {7}, {4}, {30}, {16}, {23}};
auto v3 = qvec | view::to_char;
std::cout << v3 << '\n'; // [!,(,&,$,(,%,?,1,8]

std::vector<dna4q> qcvec{{dna4::C, phred42{0}}, {dna4::A, phred42{7}}, {dna4::G, phred42{5}}, {dna4::T, phred42{3}},
                         {dna4::G, phred42{7}}, {dna4::A, phred42{4}}, {dna4::C, phred42{30}}, {dna4::T, phred42{16}}, {dna4::A, phred42{23}}};
auto v4 = qcvec | view::to_char;
std::cout << v4 << '\n'; // [C,A,G,T,G,A,C,T,A]
//! [to_char]
}

{
//! [to_rank]
dna4_vector vec = "ACTTTGATA"_dna4;
auto v = vec | view::to_rank | view::convert<unsigned>;
std::cout << v << '\n'; // [0,1,3,3,3,2,0,3,0]

std::vector<phred42> qvec{{0}, {7}, {5}, {3}, {7}, {4}, {30}, {16}, {23}};
auto v3 = qvec | view::to_rank | view::convert<unsigned>;
std::cout << v3 << '\n'; // [0,7,5,3,7,4,30,16,23]

std::vector<dna4q> qcvec{{dna4::C, phred42{0}}, {dna4::A, phred42{7}}, {dna4::G, phred42{5}}, {dna4::T, phred42{3}},
                         {dna4::G, phred42{7}}, {dna4::A, phred42{4}}, {dna4::C, phred42{30}}, {dna4::T, phred42{16}}, {dna4::A, phred42{23}}};
auto v4 = qcvec | view::to_rank | view::convert<unsigned>;
std::cout << v4 << '\n'; // [1,28,22,15,30,16,121,67,92]
//! [to_rank]
}

}
