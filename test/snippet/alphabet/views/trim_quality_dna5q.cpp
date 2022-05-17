#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/alphabet/views/trim_quality.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    std::vector<seqan3::dna5q> vec{{'A'_dna5, 'I'_phred42},
                                   {'G'_dna5, 'I'_phred42},
                                   {'G'_dna5, '?'_phred42},
                                   {'A'_dna5, '5'_phred42},
                                   {'T'_dna5, '+'_phred42}};

    // trim by phred_value
    auto v1 = vec | seqan3::views::trim_quality(20u);
    seqan3::debug_stream << v1 << '\n'; // AGGA

    // trim by quality character; in this case the nucleotide part of the character is irrelevant
    auto v2 = vec | seqan3::views::trim_quality(seqan3::dna5q{'C'_dna5, '5'_phred42});
    seqan3::debug_stream << v2 << '\n'; // AGGA

    // combinability
    auto v3 = seqan3::views::trim_quality(vec, 20u) | seqan3::views::to_char;
    seqan3::debug_stream << v3 << '\n'; // AGGA
}
