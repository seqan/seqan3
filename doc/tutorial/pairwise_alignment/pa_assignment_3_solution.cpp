#include <utility>
#include <vector>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/debug_stream.hpp>

using seqan3::operator""_aa27;

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
    auto config = seqan3::align_cfg::mode{seqan3::global_alignment} |
                  seqan3::align_cfg::scoring{seqan3::aminoacid_scoring_scheme{
                      seqan3::aminoacid_similarity_matrix::BLOSUM62}} |
                  seqan3::align_cfg::aligned_ends{seqan3::free_ends_second};

    for (auto const & res : seqan3::align_pairwise(source, config))
        seqan3::debug_stream << "Score: " << res.score() << '\n';
}
