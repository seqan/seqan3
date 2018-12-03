#include <seqan3/std/view/subrange.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

using namespace seqan3;

int main()
{
//! [example]
dna4_vector s{"ACTTTGATAN"_dna4};
auto v1 = view::subrange<decltype(begin(s)), decltype(end(s))>{begin(s) + 2, end(s)} | view::to_char; // == "TTTGATAA"
//! [example]

debug_stream << v1 << '\n';
}
