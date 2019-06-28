#include <seqan3/alphabet/composite/alphabet_variant.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>

int main()
{
    using seqan3::operator""_dna5;
    using seqan3::operator""_rna5;

    seqan3::alphabet_variant<seqan3::dna5, seqan3::gap> letter{};          // implicitly 'A'_dna5
    seqan3::alphabet_variant<seqan3::dna5, seqan3::gap> letter2{'C'_dna5}; // constructed from alternative (== 'C'_dna5)
    seqan3::alphabet_variant<seqan3::dna5, seqan3::gap> letter3{'U'_rna5}; // constructed from type that alternative is constructible from (== 'T'_dna5)

    letter2.assign_char('T');                       // == 'T'_dna5
    letter2.assign_char('-');                       // == gap{}
    letter2.assign_char('K');                       // unknown characters map to the default/unknown
                                                    // character of the first alternative type (== 'N'_dna5)

    letter2 = seqan3::gap{};                        // assigned from alternative (== gap{})
    letter2 = 'U'_rna5;                             // assigned from type that alternative is assignable from (== 'T'_dna5)

    seqan3::dna5 letter4 = letter2.convert_to<seqan3::dna5>();
}
