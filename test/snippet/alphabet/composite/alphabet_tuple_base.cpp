#include <seqan3/alphabet/composite/alphabet_tuple_base.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{

{
//! [value_construction]
qualified<dna4, phred42> letter1{'C'_dna4};    // creates {'C'_dna4, phred42{0}}
qualified<dna4, phred42> letter2{phred42{1}}; // creates {'A'_dna4, phred42{1}}
//! [value_construction]
}

{
//! [subtype_construction]
// The following creates {'C'_dna4, phre42{0}}
qualified<dna4, phred42> letter1{'C'_dna4};
// The following also creates {'C'_dna4, phred42{0}}, since rna4 assignable to dna4
qualified<dna4, phred42> letter2{'C'_rna4};

if (letter1 == letter2)
    debug_stream << "yeah\n"; // yeah
//! [subtype_construction]
}

{
//! [value_assignment]
qualified<dna4, phred42> letter1{'T'_dna4, phred42{1}};

letter1 = 'C'_dna4; // yields {'C'_dna4, phred42{1}}
letter1 = phred42{2}; // yields {'C'_dna4, phred42{2}}
//! [value_assignment]
}

{
//! [subtype_assignment]
qualified<dna4, phred42> letter1{'T'_dna4, phred42{1}};

letter1 = 'C'_rna4; // yields {'C'_dna4, phred42{1}}
//! [subtype_assignment]
}

}
