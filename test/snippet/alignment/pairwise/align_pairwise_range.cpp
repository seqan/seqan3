#include <vector>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/ranges>

int main()
{
    using seqan3::operator""_dna4;

    std::vector data1{"AGTGCTACG"_dna4, "AGTAGACTACG"_dna4, "AGTTACGAC"_dna4};
    std::vector data2{"ACGTGCGACTAG"_dna4, "ACGTACGACACG"_dna4, "AGTAGCGATCG"_dna4};

    // Configure the alignment kernel.
    auto config = seqan3::align_cfg::mode{seqan3::global_alignment} |
                  seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{}};

    // Compute the alignment over a range of pairs.
    for (auto const & res : seqan3::align_pairwise(std::view::zip(data1, data2), config))
        seqan3::debug_stream << "The score: " << res.score() << "\n";
}
