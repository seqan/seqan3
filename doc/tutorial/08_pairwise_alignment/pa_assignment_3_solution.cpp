// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <utility>
#include <vector>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    std::vector vec{"MANLGYZW"_aa27, "LCKRLGNM"_aa27, "KPSKPRDYEDG"_aa27, "EQMCITQYR"_aa27};

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
    auto config = seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{false},
                                                   seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                                                   seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
                                                   seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}}
                | seqan3::align_cfg::scoring_scheme{seqan3::aminoacid_scoring_scheme{
                    seqan3::aminoacid_similarity_matrix::blosum62}}
                | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10},
                                                     seqan3::align_cfg::extension_score{-1}};

    for (auto const & res : seqan3::align_pairwise(source, config))
        seqan3::debug_stream << "Score: " << res.score() << '\n';
}
