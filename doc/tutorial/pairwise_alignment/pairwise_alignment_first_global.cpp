#include <utility>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

using seqan3::operator""_dna4;

int main()
{
    seqan3::dna4_vector s1 = "ACGTGAACTGACT"_dna4;
    seqan3::dna4_vector s2 = "ACGAAGACCGAT"_dna4;

    // Configure the alignment kernel.
    auto config = seqan3::align_cfg::mode{seqan3::global_alignment} |
                  seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{}};

    // Invoke the pairwise alignment which returns a lazy range over alignment results.
    auto results = seqan3::align_pairwise(std::tie(s1, s2), config);
    auto & res = *results.begin();
    seqan3::debug_stream << "Score: " << res.score() << '\n';
}
