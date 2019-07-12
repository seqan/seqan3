// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <vector>

#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_mode.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "alignment_fixture.hpp"

namespace seqan3::test::alignment::fixture::semi_global::affine::banded
{

inline constexpr auto align_config = align_cfg::mode{global_alignment} |
                                     align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
                                     align_cfg::band{static_band{lower_bound{-4}, upper_bound{8}}};

inline constexpr auto align_config_semi_seq1 = align_config | align_cfg::aligned_ends{free_ends_first};
inline constexpr auto align_config_semi_seq2 = align_config | align_cfg::aligned_ends{free_ends_second};

static auto dna4_01_semi_first = []()
{
    return alignment_fixture
    {
        "TTTTTACGTATGTCCCCC"_dna4,
        "ACGTAAAACGT"_dna4,
        align_config_semi_seq1 | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
        10,
        "ACGT---ATGT",
        "ACGTAAAACGT",
        alignment_coordinate{detail::column_index_type{5u}, detail::row_index_type{0u}},
        alignment_coordinate{detail::column_index_type{13u}, detail::row_index_type{11u}}
    };
}();

// TODO: Fix the computation of the trailing gap, which depends on whether the score in the last cell
// was a gap or not (the gap must be extended in this case.) -> or throw invalid_alignment_configuration if
// the score cannot be computed.
// static auto dna4_02_semi_first = []()
// {
//     return alignment_fixture
//     {
//         "ACGTAAAACGT"_dna4,
//         "TTTTTACGTATGTCCCCC"_dna4,
//         align_config_semi_seq1 | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
//         -13,
//         "-----ACGTA--------AAACGT",
//         "TTTTTACGTATGTCCCCC------",
//         alignment_coordinate{detail::column_index_type{0u}, detail::row_index_type{0u}},
//         alignment_coordinate{detail::column_index_type{5u}, detail::row_index_type{18u}}
//     };
// }();

static auto dna4_03_semi_second = []()
{
    return alignment_fixture
    {
        "TTTTTACGTATGTCCCCC"_dna4,
        "ACGTAAAACGT"_dna4,
        align_config_semi_seq2 | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
        -19,
        "TTTTTACGTATGTCCCCC",
        "GTAAAACGT---------",
        alignment_coordinate{detail::column_index_type{0u}, detail::row_index_type{2u}},
        alignment_coordinate{detail::column_index_type{18u}, detail::row_index_type{11u}}
    };
}();

static auto dna4_04_semi_second = []()
{
    return alignment_fixture
    {
        "ACGTAAAACGT"_dna4,
        "TTTTTACGTATGTCCCCC"_dna4,
        align_config_semi_seq2 | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
        -5,
        "ACGTAAAACGT",
        "------TACGT",
        alignment_coordinate{detail::column_index_type{0u}, detail::row_index_type{4u}},
        alignment_coordinate{detail::column_index_type{11u}, detail::row_index_type{9u}}
    };
}();

} // namespace seqan3::test::alignment::fixture::global::affine::banded
