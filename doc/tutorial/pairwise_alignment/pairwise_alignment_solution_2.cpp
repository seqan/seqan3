#include <utility>
#include <vector>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/pairwise_combine.hpp>

using seqan3::operator""_dna4;

int main()
{
    std::vector vec{"ACGTGAACTGACT"_dna4,
                    "ACGAAGACCGAT"_dna4,
                    "ACGTGACTGACT"_dna4,
                    "AGGTACGAGCGACACT"_dna4};

    // Configure the alignment kernel.
    auto config = seqan3::align_cfg::method_global{
                      seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                      seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                      seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                      seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}} |
                  seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}};

    for (auto const & res : seqan3::align_pairwise(seqan3::views::pairwise_combine(vec), config))
        seqan3::debug_stream << "Score: " << res.score() << '\n';
}
