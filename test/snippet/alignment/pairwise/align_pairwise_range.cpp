#include <vector>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

int main()
{
    std::vector data1{"AGTGCTACG"_dna4, "AGTAGACTACG"_dna4, "AGTTACGAC"_dna4};
    std::vector data2{"ACGTGCGACTAG"_dna4, "ACGTACGACACG"_dna4, "AGTAGCGATCG"_dna4};

    // Configure the alignment kernel.
    auto config = align_cfg::mode{global_alignment} |
                  align_cfg::scoring{nucleotide_scoring_scheme{}};

    // Compute the alignment over a range of pairs.
    for (auto const & res : align_pairwise(std::view::zip(data1, data2), config))
        debug_stream << "The score: " << res.score() << "\n";
}
