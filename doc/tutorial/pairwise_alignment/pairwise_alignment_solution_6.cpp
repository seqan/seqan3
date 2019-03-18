#include <utility>
#include <vector>

#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/range/view/pairwise_combine.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

int main()
{
    std::vector vec{"ACGTGACTGACT"_dna4,
                    "ACGAAGACCGAT"_dna4,
                    "ACGTGACTGACT"_dna4,
                    "AGGTACGAGCGACACT"_dna4};

    using pair_t = decltype(std::tie(vec[0], vec[0]));
    std::vector<pair_t> source;
    for (unsigned i = 0; i < vec.size() - 1; ++i)
    {
        for (unsigned j = i + 1; j < vec.size(); ++j)
        {
            source.push_back(std::tie(vec[i], vec[j]));
        }
    }

    // Configure the alignment kernel.
    auto config = align_cfg::edit |
                  align_cfg::max_error{7u} |
                  align_cfg::result{with_alignment};

    auto filter_v = std::view::filter([](auto && res) { return res.score() > -7;});

    for (auto const & res : align_pairwise(view::pairwise_combine(vec), config) | view::persist | filter_v)
    {
        debug_stream << "Score: " << res.score() << '\n';
        debug_stream << res.alignment() << '\n';
    }
}
