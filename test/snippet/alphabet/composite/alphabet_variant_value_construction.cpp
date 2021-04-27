#include <seqan3/alphabet/composite/alphabet_variant.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::alphabet_variant<seqan3::dna4, seqan3::gap> letter1{'C'_dna4}; // or
    seqan3::alphabet_variant<seqan3::dna4, seqan3::gap> letter2 = seqan3::gap{};
}
