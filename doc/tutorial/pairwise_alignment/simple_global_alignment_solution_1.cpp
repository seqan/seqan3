#include <utility>
#include <vector>

#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

using namespace seqan3;

int main()
{
    std::vector vec{"ACGTGAACTGACT"_dna4,
                    "ACGAAGACCGAT"_dna4,
                    "ACGTGACTGACT"_dna4,
                    "AGGTACGAGCGACACT"_dna4};

    using pair_t = decltype(std::tie(vec[0], vec[0]));
    std::vector<pair_t> source;
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        for (unsigned j = i + 1; j < vec.size(); ++j)
        {
            source.push_back(std::tie(vec[i], vec[j]));
        }
    }

    // Configure the alignment kernel.
    auto config = align_cfg::mode{global_alignment} |
                  align_cfg::scoring{nucleotide_scoring_scheme{}};

    for (auto const & res : align_pairwise(source, config))
        debug_stream << "Score: " << res.score() << '\n';
}
