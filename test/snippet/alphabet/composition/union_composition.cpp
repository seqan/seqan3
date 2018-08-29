#include <seqan3/alphabet/composition/union_composition.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <gtest/gtest.h>

using namespace seqan3;
int main()
{

{
//! [usage]
union_composition<dna5, gap> letter{};         // implicitly dna5::A
union_composition<dna5, gap> letter2{dna5::C}; // constructed from alternative (== dna5::C)
union_composition<dna5, gap> letter3{rna5::U}; // constructed from type that alternative is constructable from (== dna5::T)

letter2.assign_char('T');                      // == dna5::T
letter2.assign_char('-');                      // == gap::GAP
letter2.assign_char('K');                      // unknown characters map to the default/unknown
                                               // character of the first alternative type (== dna5::N)

letter2 = gap::GAP;                            // assigned from alternative (== gap::GAP)
letter2 = rna5::U;                             // assigned from type that alternative is assignable from (== dna5::T)

dna5 letter4 = letter2.convert_to<dna5>();     // this works
// gap letter5  = letter2.convert_to<gap>();   // this throws an exception, because the set value was dna5::T
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
union_composition<dna4, gap> letter1{dna4::C}; // or
union_composition<dna4, gap> letter2 = gap::GAP;
//! [value construction]
(void) letter1;
(void) letter2;
}

{
//! [conversion]
union_composition<dna4, gap> letter1{rna4::C};
//! [conversion]
(void) letter1;
}

{
//! [subtype_construction]
union_composition<dna4, gap> letter1{};
letter1 = rna4::C;
//! [subtype_construction]
}

}
