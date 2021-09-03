// generated from test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@_inherit.cpp.in

//![main]
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

struct my_rna4 : public seqan3::rna4
{
    // using seqan3::rna4::rna4; // uncomment to import implicit conversion shown by letter1
};

struct my_dna4 : public seqan3::dna4
{};

int main()
{
    using namespace seqan3::literals;

    // my_rna4 letter1 = 'C'_dna4; // NO automatic implicit conversion!
    // seqan3::rna4 letter2 = my_dna4{}; // seqan3::rna4 only allows implicit conversion from seqan3::dna4!
}
//![main]

#include <seqan3/utility/concept/exposition_only/core_language.hpp>

static_assert(seqan3::implicitly_convertible_to<seqan3::dna4, seqan3::rna4>);
static_assert(!seqan3::implicitly_convertible_to<seqan3::dna4, my_rna4>);
static_assert(!seqan3::implicitly_convertible_to<my_dna4, seqan3::rna4>);
