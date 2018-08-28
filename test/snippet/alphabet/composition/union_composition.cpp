#include <seqan3/alphabet/composition/union_composition.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <gtest/gtest.h>

using namespace seqan3;
int main()
{

{
//! [variant]
union_composition<dna4, gap> my_letter{};
union_composition<dna4, gap> converted_letter{dna4::C};
// doesn't work:
// union_composition<dna4, gap> my_letter{'A'};
union_composition<dna4, gap>{}.assign_char('C'); // <- this does!
union_composition<dna4, gap>{}.assign_char('-'); // gap character
union_composition<dna4, gap>{}.assign_char('K'); // unknown characters map to the default/unknown
                                                 // character of the first alphabet type (i.e. A of dna4)
if (my_letter.to_char() == 'A')
std::cout << "yeah\n"; // "yeah";
//! [variant]
}

{
//! [construct base]
using alphabet_t = union_composition<dna4, dna5, gap>;

constexpr alphabet_t letter0{gap::GAP};
constexpr alphabet_t letter1 = dna5::C;
constexpr alphabet_t letter2 = {dna4::G};
constexpr alphabet_t letter3 = static_cast<alphabet_t>(dna4::T);

assert(letter0.to_rank() == 9);
assert(letter1.to_rank() == 5);
assert(letter2.to_rank() == 2);
assert(letter3.to_rank() == 3);
//! [construct base]
}

{
//! [assign base]
using alphabet_t = union_composition<dna4, dna5, gap>;

alphabet_t letter;

letter = dna5::A;
assert(letter.to_rank() == 4);

letter = {dna5::C};
assert(letter.to_rank() == 5);

letter = static_cast<alphabet_t>(dna5::G);
assert(letter.to_rank() == 6);
//! [assign base]
}

{
//! [has_alternative]
using union_t = union_composition<dna4, gap>;

static_assert(union_t::has_alternative<dna4>(), "should be true");
static_assert(union_t::has_alternative<gap>(), "should be true");
static_assert(!union_t::has_alternative<dna5>(), "should be false");
//! [has_alternative]
}

{
//! [value construction]
union_composition<dna4, gap> letter1{dna4::C}; // or
union_composition<dna4, gap> letter2 = gap::GAP;
//! [value construction]
(void) letter2;
}

{
//! [reoccurring construction]
using alphabet_t = union_composition<dna4, dna4>;

constexpr alphabet_t letter0{std::in_place_index_t<0>{}, dna4::A};
constexpr alphabet_t letter4{std::in_place_index_t<1>{}, dna4::A};

EXPECT_EQ(letter0.to_rank(), 0);
EXPECT_EQ(letter4.to_rank(), 4);
//! [reoccurring construction]
}

{
//! [assign by base]
union_composition<dna4, gap> letter1{};
letter1 = rna4::A;
//! [assign by base]
}

{
//! [conversion]
union_composition<dna4, gap> letter1{rna4::C};
//! [conversion]
}

{
//! [subtype_construction]
union_composition<dna4, gap> letter1{};
letter1 = rna4::C;
//! [subtype_construction]
}

}
