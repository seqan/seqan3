#include <vector>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/ranges>

int main()
{
    using namespace seqan3::literals;

    std::vector data1{"AGTGCTACG"_dna4, "AGTAGACTACG"_dna4, "AGTTACGAC"_dna4};
    std::vector data2{"ACGTGCGACTAG"_dna4, "ACGTACGACACG"_dna4, "AGTAGCGATCG"_dna4};

    // Configure the alignment kernel.
    auto config = seqan3::align_cfg::method_global{} |
                  seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}};

    // Compute the alignment over a range of pairs.
    for (auto const & res : seqan3::align_pairwise(seqan3::views::zip(data1, data2), config))
        seqan3::debug_stream << "The score: " << res.score() << "\n";
}
