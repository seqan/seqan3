#include <utility>
#include <vector>

#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
    std::vector vec{"MANLGYZW"_aa27,
                    "LCKRLGNM"_aa27,
                    "KPSKPRDYEDG"_aa27,
                    "EQMCITQYR"_aa27};

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
                  align_cfg::scoring{aminoacid_scoring_scheme{aminoacid_similarity_matrix::BLOSUM62}} |
                  align_cfg::aligned_ends{free_ends_second};

    for (auto const & res : align_pairwise(source, config))
        debug_stream << "Score: " << res.score() << '\n';
}
