// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

//! [include]
#include <seqan3/alignment/configuration/all.hpp>
//! [include]

//! [include_scoring_scheme]
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/hamming_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
//! [include_scoring_scheme]

//! [include_gap_cost_affine]
#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
//! [include_gap_cost_affine]

//! [include_method]
#include <seqan3/alignment/configuration/align_config_method.hpp>
//! [include_method]

//! [include_output]
#include <seqan3/alignment/configuration/align_config_output.hpp>
//! [include_output]

//! [include_band]
#include <seqan3/alignment/configuration/align_config_band.hpp>
//! [include_band]

//! [include_edit]
#include <seqan3/alignment/configuration/align_config_edit.hpp>
//! [include_edit]

#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>

int main()
{
    {
        //! [method_global_free_end_gaps]
        // Example of a semi-global alignment where leading and trailing gaps in the
        // second sequence are not penalised:
        auto config = seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{false},
                                                       seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                                                       seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
                                                       seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}};
        //! [method_global_free_end_gaps]
    }

    {
        //! [scoring_scheme]
        using namespace seqan3::literals;

        // Define a simple scoring scheme with match and mismatch cost and get the score.
        seqan3::nucleotide_scoring_scheme nc_scheme{seqan3::match_score{4}, seqan3::mismatch_score{-5}};
        auto sc_nc = nc_scheme.score('A'_dna4, 'C'_dna4); // sc_nc == -5.

        // Define a amino acid similarity matrix and get the score.
        seqan3::aminoacid_scoring_scheme aa_scheme{};
        aa_scheme.set_similarity_matrix(seqan3::aminoacid_similarity_matrix::blosum30);
        auto sc_aa = aa_scheme.score('M'_aa27, 'K'_aa27); // sc_aa == 2.
        //! [scoring_scheme]
    }

    {
        //! [gap_cost_affine]

        // Define a gap scheme with custom gap scores.
        seqan3::align_cfg::gap_cost_affine affine_scheme{seqan3::align_cfg::open_score{-10},
                                                         seqan3::align_cfg::extension_score{-1}};

        int open_score = affine_scheme.open_score;           // == -10
        int extension_score = affine_scheme.extension_score; // == -1
        //! [gap_cost_affine]
    }

    {
        //! [output]

        // Configure the alignment to only compute the score.
        auto cfg = seqan3::align_cfg::output_score{};
        //! [output]
    }

    {
        //! [band]

        // Configure a banded alignment.
        auto cfg = seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-4},
                                                      seqan3::align_cfg::upper_diagonal{4}};
        //! [band]
    }

    {
        //! [edit]

        // Configure an edit distance alignment.
        auto cfg = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme;
        //! [edit]
    }
}
