#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/range/view/trim.hpp>
#include <seqan3/range/view/to_char.hpp>

int main()
{
    using seqan3::operator""_dna5;
    using seqan3::operator""_phred42;
    std::vector<seqan3::dna5q> vec{{'A'_dna5, 'I'_phred42},
                                   {'G'_dna5, 'I'_phred42},
                                   {'G'_dna5, '?'_phred42},
                                   {'A'_dna5, '5'_phred42},
                                   {'T'_dna5, '+'_phred42}};
    std::vector<seqan3::dna5q> cmp{{'A'_dna5, 'I'_phred42},
                                   {'G'_dna5, 'I'_phred42},
                                   {'G'_dna5, '?'_phred42},
                                   {'A'_dna5, '5'_phred42}};

    // trim by phred_value
    auto v1 = vec | seqan3::view::trim(20u);
    assert(std::vector<seqan3::dna5q>(v1) == cmp);

    // trim by quality character; in this case the nucleotide part of the character is irrelevant
    auto v2 = vec | seqan3::view::trim(seqan3::dna5q{'C'_dna5, '5'_phred42});
    assert(std::vector<seqan3::dna5q>(v2) == cmp);

    // combinability
    auto v3 = seqan3::view::trim(vec, 20u) | seqan3::view::to_char;
    assert(std::ranges::equal(std::string_view{"AGGA"}, v3));
}
