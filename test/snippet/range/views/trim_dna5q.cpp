#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/range/views/trim.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/to.hpp>

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
    auto v1 = vec | seqan3::views::trim(20u);
    assert((v1 | seqan3::views::to<std::vector>) == cmp);

    // trim by quality character; in this case the nucleotide part of the character is irrelevant
    auto v2 = vec | seqan3::views::trim(seqan3::dna5q{'C'_dna5, '5'_phred42});
    assert((v2 | seqan3::views::to<std::vector>) == cmp);

    // combinability
    auto v3 = seqan3::views::trim(vec, 20u) | seqan3::views::to_char;
    assert(std::ranges::equal(std::string{"AGGA"}, v3 | seqan3::views::to<std::string>));
}
