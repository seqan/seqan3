#include <seqan3/alphabet/composite/alphabet_variant.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <gtest/gtest.h>

using namespace seqan3;
int main()
{

{
//! [usage]
alphabet_variant<dna5, gap> letter{};          // implicitly 'A'_dna5
alphabet_variant<dna5, gap> letter2{'C'_dna5}; // constructed from alternative (== 'C'_dna5)
alphabet_variant<dna5, gap> letter3{'U'_rna5}; // constructed from type that alternative is constructable from (== 'T'_dna5)

letter2.assign_char('T');                       // == 'T'_dna5
letter2.assign_char('-');                       // == gap{}
letter2.assign_char('K');                       // unknown characters map to the default/unknown
                                                // character of the first alternative type (== 'N'_dna5)

letter2 = gap{};                                // assigned from alternative (== gap{})
letter2 = 'U'_rna5;                             // assigned from type that alternative is assignable from (== 'T'_dna5)

dna5 letter4 = letter2.convert_to<dna5>();      // this works
// gap letter5  = letter2.convert_to<gap>();    // this throws an exception, because the set value was 'T'_dna5
//! [usage]
(void) letter4;
}

{
//! [char_representation]
alphabet_variant<dna4, dna5> var;
var.assign_char('A');             // will be in the "dna4-state"
var = 'A'_dna5;                   // will be in the "dna5-state"
//! [char_representation]
}

{
//! [holds_alternative]
using variant_t = alphabet_variant<dna5, gap>;

static_assert(variant_t::holds_alternative<dna5>(), "dna5 is an alternative of variant_t");
static_assert(!variant_t::holds_alternative<dna4>(), "dna4 is not an alternative of variant_t");
static_assert(variant_t::holds_alternative<gap>(), "gap is an alternative of variant_t");
//! [holds_alternative]
}

{
//! [value construction]
alphabet_variant<dna4, gap> letter1{'C'_dna4}; // or
alphabet_variant<dna4, gap> letter2 = gap{};
//! [value construction]
(void) letter1;
(void) letter2;
}

{
//! [conversion]
alphabet_variant<dna4, gap> letter1{'C'_rna4};
//! [conversion]
(void) letter1;
}

{
//! [subtype_construction]
alphabet_variant<dna4, gap> letter1{};
letter1 = 'C'_rna4;
//! [subtype_construction]
}

}
