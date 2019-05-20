#include <utility>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
    auto seq1 = "TTACGTACGGACTAGCTACAACATTACGGACTAC"_dna4;
    auto seq2 = "GGACGACATGACGTACGACTTTACGTACGACTAGC"_dna4;

    // Configure the alignment kernel.
    auto config = align_cfg::mode{global_alignment} |
                  align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-2}}} |
                  align_cfg::gap{gap_scheme{gap_score{-4}}} |
                  align_cfg::aligned_ends{free_ends_all} |
                  align_cfg::result{with_alignment};

    for (auto const & res : align_pairwise(std::tie(seq1, seq2), config))
    {
        debug_stream << "Score: " << res.score() << '\n';
        debug_stream << "Begin: " << res.front_coordinate() << '\n';
        debug_stream << "End: " << res.back_coordinate() << '\n';
        debug_stream << "Alignment: \n" << res.alignment() << '\n';
    }
}
