#include <seqan3/range/view/get.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>

using namespace seqan3;

int main()
{
//! [usage]
// Create a vector of dna4 quality composition alphabet.
std::vector<dna4q> qv{{dna4::A, phred42{0}}, {dna4::C, phred42{1}}, {dna4::G, phred42{2}}, {dna4::T, phred42{3}}};

std::cout << (qv | view::get<0> | view::to_char) << std::endl; // Prints [A,C,G,T]
std::cout << (qv | view::get<1> | view::to_char) << std::endl; // Prints [!,",#,$]
//! [usage]
}
