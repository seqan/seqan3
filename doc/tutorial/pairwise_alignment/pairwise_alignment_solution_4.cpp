#include <utility>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

using seqan3::operator""_dna4;

int main()
{
    auto seq1 = "TTACGTACGGACTAGCTACAACATTACGGACTAC"_dna4;
    auto seq2 = "GGACGACATGACGTACGACTTTACGTACGACTAGC"_dna4;

    // Configure the alignment kernel.
    auto config = seqan3::align_cfg::mode{seqan3::global_alignment} |
                  seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{
                      seqan3::match_score{4}, seqan3::mismatch_score{-2}}} |
                  seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-4}}} |
                  seqan3::align_cfg::aligned_ends{seqan3::free_ends_all} |
                  seqan3::align_cfg::result{seqan3::with_alignment};

    for (auto const & res : seqan3::align_pairwise(std::tie(seq1, seq2), config))
    {
        seqan3::debug_stream << "Score: " << res.score() << '\n';
        seqan3::debug_stream << "Begin: " << res.front_coordinate() << '\n';
        seqan3::debug_stream << "End: " << res.back_coordinate() << '\n';
        seqan3::debug_stream << "Alignment: \n" << res.alignment() << '\n';
    }
}
