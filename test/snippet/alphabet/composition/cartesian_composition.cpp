#include <seqan3/alphabet/composition/cartesian_composition.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>

using namespace seqan3;

int main()
{
{
//! [value_construction]
cartesian_composition<dna4, aa27> letter1{dna4::C}; // creates {dna4::C, aa27::A}
cartesian_composition<dna4, aa27> letter2{aa27::K}; // creates {dna4::A, aa27::K}
//! [value_construction]
}

{
//! [subtype_construction]
// The following creates {dna4::C, aa27::A}
cartesian_composition<gapped<dna4>, aa27> letter1{dna4::C};
// The following creates {dna4::A, dna4::C}, since dna4 is also a type in the list
cartesian_composition<gapped<dna4>, dna4> letter2{dna4::C};
// The following creates {dna5::C, dna15::A}, since dna5 is the first type assignable to dna4
cartesian_composition<dna5, dna15> letter2{dna4::C};
//! [subtype_construction]
}

{
//! [value_assignment]
cartesian_composition<dna4, aa27> letter1{dna4::T, aa27::K};

letter1 = dna4::C // yields {dna4::C, aa27::K}
letter1 = aa27::F // yields {dna4::C, aa27::F}
//! [value_assignment]
}

{
//! [subtype_assignment]
cartesian_composition<dna4, aa27> letter1{dna4::T, aa27::K};

letter1 = rna4::C; // yields {dna4::C, aa27::K}
//! [subtype_assignment]
}
}
