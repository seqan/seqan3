#include <seqan3/alphabet/composition/cartesian_composition.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

using namespace seqan3;

int main()
{

{
//! [value_construction]
qualified<dna4, phred42> letter1{dna4::C};    // creates {dna4::C, phred42{0}}
qualified<dna4, phred42> letter2{phred42{1}}; // creates {dna4::A, phred42{1}}
//! [value_construction]
}

{
//! [subtype_construction]
// The following creates {dna4::C, phre42{0}}
qualified<dna4, phred42> letter1{dna4::C};
// The following also creates {dna4::C, phred42{0}}, since rna4 assignable to dna4
qualified<dna4, phred42> letter2{rna4::C};

if (letter1 == letter2)
    debug_stream << "yeah\n"; // yeah
//! [subtype_construction]
}

{
//! [value_assignment]
qualified<dna4, phred42> letter1{dna4::T, phred42{1}};

letter1 = dna4::C; // yields {dna4::C, phred42{1}}
letter1 = phred42{2}; // yields {dna4::C, phred42{2}}
//! [value_assignment]
}

{
//! [subtype_assignment]
qualified<dna4, phred42> letter1{dna4::T, phred42{1}};

letter1 = rna4::C; // yields {dna4::C, phred42{1}}
//! [subtype_assignment]
}

}
