#include <utility>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
    auto seq1 = "QFSEEILSDIYCWMLQCGQERAV"_aa27;
    auto seq2 = "AFLPGWQEENKLSKIWMKDCGCLW"_aa27;

    // Configure the alignment kernel.
    auto config = align_cfg::mode{global_alignment} |
                  align_cfg::scoring{aminoacid_scoring_scheme{aminoacid_similarity_matrix::BLOSUM62}} |
                  align_cfg::gap{gap_scheme{gap_score{-2}, gap_open_score{-9}}};

    for (auto const & res : align_pairwise(std::tie(seq1, seq2), config))
        debug_stream << "Score: " << res.score() << '\n';
}
