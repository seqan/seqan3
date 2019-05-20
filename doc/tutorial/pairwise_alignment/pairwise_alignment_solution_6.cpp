#include <utility>
#include <vector>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/pairwise_combine.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

int main()
{
    std::vector vec{"ACGTGACTGACT"_dna4,
                    "ACGAAGACCGAT"_dna4,
                    "ACGTGACTGACT"_dna4,
                    "AGGTACGAGCGACACT"_dna4};

    // Configure the alignment kernel.
    auto config = align_cfg::edit |
                  align_cfg::max_error{7u} |
                  align_cfg::result{with_score};

    auto filter_v = std::view::filter([](auto && res) { return res.score() >= -6;});

    for (auto const & res : align_pairwise(view::pairwise_combine(vec), config) | view::persist | filter_v)
    {
        debug_stream << "Score: " << res.score() << '\n';
    }
}
