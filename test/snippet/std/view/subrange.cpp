#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

int main()
{
//! [example]
dna4_vector s{"ACTTTGATAA"_dna4};
using iterator = dna4_vector::iterator;
auto v1 = std::ranges::subrange<iterator, iterator>{begin(s) + 2, end(s)} | view::to_char; // == "TTTGATAA"
//! [example]

debug_stream << v1 << '\n';
}
