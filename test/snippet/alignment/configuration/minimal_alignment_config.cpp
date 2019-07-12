#include <seqan3/alignment/configuration/align_config_mode.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

using seqan3::operator""_dna4;

int main()
{
    auto min_cfg = seqan3::align_cfg::mode{seqan3::global_alignment} |
                   seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4}, seqan3::mismatch_score{-5}}};

    auto seq1 = "ACGT"_dna4;
    auto seq2 = "ACCT"_dna4;
    for (auto res : align_pairwise(std::tie(seq1, seq2), min_cfg))
        seqan3::debug_stream << res.score() << '\n';
}
