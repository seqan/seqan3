#include <seqan3/alignment/configuration/align_config_mode.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

using namespace seqan3;

int main()
{
    auto min_cfg = align_cfg::mode{global_alignment} |
                   align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}};

    auto seq1 = "ACGT"_dna4;
    auto seq2 = "ACCT"_dna4;
    for (auto res : align_pairwise(std::tie(seq1, seq2), min_cfg))
        debug_stream << res.score() << '\n';
}
