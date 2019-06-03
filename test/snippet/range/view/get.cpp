#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/get.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>

using namespace seqan3;

int main()
{
//! [usage]
// Create a vector of dna4 quality composite alphabet.
std::vector<dna4q> qv{{'A'_dna4, phred42{0}}, {'C'_dna4, phred42{1}}, {'G'_dna4, phred42{2}}, {'T'_dna4, phred42{3}}};

debug_stream << (qv | view::get<0> | view::to_char) << std::endl; // Prints [A,C,G,T]
debug_stream << (qv | view::get<1> | view::to_char) << std::endl; // Prints [!,",#,$]
//! [usage]
}
