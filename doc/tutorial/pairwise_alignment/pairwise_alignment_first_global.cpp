#include <utility>

#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
    dna4_vector s1 = "ACGTGAACTGACT"_dna4;
    dna4_vector s2 = "ACGAAGACCGAT"_dna4;

    // Configure the alignment kernel.
    auto config = align_cfg::mode{global_alignment} |
                  align_cfg::scoring{nucleotide_scoring_scheme{}};

    // Invoke the pairwise alignment which returns a lazy range over alignment results.
    auto results = align_pairwise(std::tie(s1, s2), config);
    auto & res = *results.begin();
    debug_stream << "Score: " << res.score() << '\n';
}
