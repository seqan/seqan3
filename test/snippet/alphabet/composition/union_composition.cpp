#include <seqan3/alphabet/composition/union_composition.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <gtest/gtest.h>

using namespace seqan3;
int main()
{

{
//! [usage]
union_composition<dna5, gap> letter{};          // implicitly 'A'_dna5
union_composition<dna5, gap> letter2{'C'_dna5}; // constructed from alternative (== 'C'_dna5)
union_composition<dna5, gap> letter3{'U'_rna5}; // constructed from type that alternative is constructable from (== 'T'_dna5)

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
//! [holds_alternative]
using union_t = union_composition<dna5, gap>;

static_assert(union_t::holds_alternative<dna5>(), "dna5 is an alternative of union_t");
static_assert(!union_t::holds_alternative<dna4>(), "dna4 is not an alternative of union_t");
static_assert(union_t::holds_alternative<gap>(), "gap is an alternative of union_t");
//! [holds_alternative]
}

{
//! [value construction]
union_composition<dna4, gap> letter1{'C'_dna4}; // or
union_composition<dna4, gap> letter2 = gap{};
//! [value construction]
(void) letter1;
(void) letter2;
}

{
//! [conversion]
union_composition<dna4, gap> letter1{'C'_rna4};
//! [conversion]
(void) letter1;
}

{
//! [subtype_construction]
union_composition<dna4, gap> letter1{};
letter1 = 'C'_rna4;
//! [subtype_construction]
}

}
