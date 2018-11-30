#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/range/view/to_char.hpp>

using namespace seqan3;

int main()
{
//! [usage]
// this doesn't work:
//auto v = "ACGT"_dna4 | view::to_char;

// this works, because we explicitly create an l-value of our dna vector:
auto vec = "ACGT"_dna4;
auto v = vec | view::to_char;

// using view::persist you can bind the temporary directly:
auto v2 = "ACGT"_dna4 | view::persist | view::to_char;

// note that view::persist must follow immediately after the temporary, the following does not work:
//auto v3 = "ACGT"_dna4 | view::to_char | view::persist;

// thus the function notation might be more intuitive:
auto v3 = view::persist("ACGT"_dna4) | view::to_char;
//! [usage]

(void)v;
(void)v2;
(void)v3;
}
