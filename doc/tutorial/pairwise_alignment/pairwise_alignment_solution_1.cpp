#include <utility>
#include <vector>

#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/pairwise_combine.hpp>

using namespace seqan3;

int main()
{
    std::vector vec{"ACGTGAACTGACT"_dna4,
                    "ACGAAGACCGAT"_dna4,
                    "ACGTGACTGACT"_dna4,
                    "AGGTACGAGCGACACT"_dna4};

    // Configure the alignment kernel.
    auto config = align_cfg::mode{global_alignment} |
                  align_cfg::scoring{nucleotide_scoring_scheme{}};

    for (auto const & res : align_pairwise(view::pairwise_combine(vec), config))
        debug_stream << "Score: " << res.score() << '\n';
}
