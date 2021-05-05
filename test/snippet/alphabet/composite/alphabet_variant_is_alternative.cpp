#include <seqan3/alphabet/composite/alphabet_variant.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <gtest/gtest.h>

int main()
{
    using variant_t = seqan3::alphabet_variant<seqan3::dna5, seqan3::gap>;

    static_assert(variant_t::is_alternative<seqan3::dna5>(), "dna5 is an alternative of variant_t");
    static_assert(!variant_t::is_alternative<seqan3::dna4>(), "dna4 is not an alternative of variant_t");
    static_assert(variant_t::is_alternative<seqan3::gap>(), "gap is an alternative of variant_t");
}
