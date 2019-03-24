#include <vector>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3;

int main()
{
    std::vector vec{std::pair{"AGTGCTACG"_dna4, "ACGTGCGACTAG"_dna4},
                    std::pair{"AGTAGACTACG"_dna4, "ACGTACGACACG"_dna4},
                    std::pair{"AGTTACGAC"_dna4, "AGTAGCGATCG"_dna4}};

    // Compute the alignment of a single pair.
    for (auto const & res : align_pairwise(std::tie(vec[0].first, vec[0].second), align_cfg::edit))
        debug_stream << "The score: " << res.score() << "\n";

    // Compute the alignment over a range of pairs.
    for (auto const & res : align_pairwise(vec, align_cfg::edit))
        debug_stream << "The score: " << res.score() << "\n";
}
