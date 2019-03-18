#include <seqan3/range/view/complement.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

int main()
{
//! [usage]
dna5_vector foo{"ACGTA"_dna5};

// pipe notation
auto v = foo | view::complement;                                  // == "TGCAT"

// function notation
dna5_vector v2(view::complement(foo));                            // == "TGCAT"

// generate the reverse complement:
dna5_vector v3 = foo | view::complement | std::view::reverse;          // == "TACGT"
//! [usage]
(void) v;
(void) v2;
(void) v3;
}
