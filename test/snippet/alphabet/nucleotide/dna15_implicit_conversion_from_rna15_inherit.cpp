// generated from test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@_inherit.cpp.in

//![main]
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>

struct my_dna15 : public seqan3::dna15
{
    // using seqan3::dna15::dna15; // uncomment to import implicit conversion shown by letter1
};

struct my_rna15 : public seqan3::rna15
{};

int main()
{
    using namespace seqan3::literals;

    // my_dna15 letter1 = 'C'_rna15; // NO automatic implicit conversion!
    // seqan3::dna15 letter2 = my_rna15{}; // seqan3::dna15 only allows implicit conversion from seqan3::rna15!
}
//![main]

#include <seqan3/utility/concept/exposition_only/core_language.hpp>

static_assert(seqan3::implicitly_convertible_to<seqan3::rna15, seqan3::dna15>);
static_assert(!seqan3::implicitly_convertible_to<seqan3::rna15, my_dna15>);
static_assert(!seqan3::implicitly_convertible_to<my_rna15, seqan3::dna15>);
