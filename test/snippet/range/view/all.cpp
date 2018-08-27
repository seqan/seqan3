#include <seqan3/range/view/all.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/std/view/reverse.hpp>

using namespace seqan3;
using namespace seqan3::literal;

int main()
{

{
//! [notation]
dna4_vector vec{"ACGGTC"_dna4};

// these are synonymous:
auto vec_view1 = vec | view::complement;
auto vec_view2 = view::complement(vec);

// both views "behave" like a collection of the elements 'T', 'G', 'C', 'C', 'A', 'G'
// but can be copied cheaply et cetera
//! [notation]
(void) vec_view1;
//! [retransform]
// just re-assign to a container
dna4_vector complemented = vec_view2;
assert(complemented == "TGCCAG"_dna4);

// or immediately create on container
dna4_vector reversed = vec | view::reverse;
assert(reversed == "CTGGCA"_dna4);
//! [retransform]

//! [composability]
// views can be composed iteratively
auto vec_view3 = vec | view::reverse;
auto vec_view4 = vec_view3 | view::complement;

// or in one line similar to the unix command line
auto vec_view5 = vec | view::complement | view::reverse;

// vec_view4 and vec_view5 are the reverse complement of "ACGGTC": "GACCGT"
//! [composability]
(void) vec_view4;
(void) vec_view5;
}
}
