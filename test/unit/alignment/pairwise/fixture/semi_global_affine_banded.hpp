// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
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

inline constexpr auto align_config = align_cfg::mode{align_cfg::global_alignment} |
                                     align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
                                     align_cfg::band{static_band{lower_bound{-4}, upper_bound{8}}};

inline constexpr auto align_config_semi_seq1 = align_config | align_cfg::aligned_ends{align_cfg::seq1_ends_free};
inline constexpr auto align_config_semi_seq2 = align_config | align_cfg::aligned_ends{align_cfg::seq2_ends_free};

static auto dna4_01 = []()
{
    return alignment_fixture
    {
        "TTTTTACGTATGTCCCCC"_dna4,
        "ACGTAAAACGT"_dna4,
        align_config_semi_seq1 | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
        10,
        "TTTTTACGT---ATGTCCCCC",
        "-----ACGTAAAACGT-----",
        alignment_coordinate{detail::column_index_type{0u}, detail::row_index_type{0u}},
        alignment_coordinate{detail::column_index_type{13u}, detail::row_index_type{11u}},
        std::vector
        {
        //     e,  A,  A,  C,  C,  G,  G,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
        /*A*/  1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
        /*C*/  2,  1,  1,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
        /*G*/  3,  2,  2,  2,  2,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13,
        /*T*/  4,  3,  3,  3,  3,  3,  3,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
        /*A*/  5,  4,  3,  4,  4,  4,  4,  4,  4,  4,  5,  6,  7,  8,  9, 10, 11,
        /*C*/  6,  5,  4,  3,  4,  5,  5,  5,  5,  5,  5,  5,  6,  7,  8,  9, 10,
        /*G*/  7,  6,  5,  4,  4,  4,  5,  6,  6,  6,  6,  6,  6,  6,  7,  8,  9,
        /*T*/  8,  7,  6,  5,  5,  5,  5,  5,  6,  7,  7,  7,  7,  7,  7,  7,  8,
        /*A*/  9,  8,  7,  6,  6,  6,  6,  6,  6,  6,  7,  8,  8,  8,  8,  8,  8
        },
        std::vector
        {
        //     e,  A,  A,  C,  C,  G,  G,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/NON,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,
        /*A*/U  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,
        /*C*/U  ,U  ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL ,DL ,L  ,L  ,L  ,L  ,
        /*G*/U  ,U  ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL ,DL ,L  ,L  ,
        /*T*/U  ,U  ,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL ,DL ,
        /*A*/U  ,DU ,D  ,DUL,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,
        /*C*/U  ,U  ,U  ,D  ,DL ,DUL,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,
        /*G*/U  ,U  ,U  ,U  ,D  ,D  ,DL ,DUL,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,
        /*T*/U  ,U  ,U  ,U  ,DU ,DU ,D  ,D  ,DL ,DUL,DU ,DU ,DU ,DU ,D  ,D  ,DL ,
        /*A*/U  ,DU ,DU ,U  ,DU ,DU ,DU ,DU ,D  ,D  ,DL ,DUL,DU ,DU ,DU ,DU ,D
        }
    };
}();

static auto dna4_02 = []()
{
    return alignment_fixture
    {
        "ACGTAAAACGT"_dna4,
        "TTTTTACGTATGTCCCCC"_dna4,
        align_config_semi_seq1 | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
        -13,
        "-----ACGTA--------AAACGT",
        "TTTTTACGTATGTCCCCC------",
        alignment_coordinate{detail::column_index_type{0u}, detail::row_index_type{0u}},
        alignment_coordinate{detail::column_index_type{5u}, detail::row_index_type{18u}},
        std::vector
        {
        //     e,  A,  A,  C,  C,  G,  G,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
        /*A*/  1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
        /*C*/  2,  1,  1,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
        /*G*/  3,  2,  2,  2,  2,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13,
        /*T*/  4,  3,  3,  3,  3,  3,  3,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
        /*A*/  5,  4,  3,  4,  4,  4,  4,  4,  4,  4,  5,  6,  7,  8,  9, 10, 11,
        /*C*/  6,  5,  4,  3,  4,  5,  5,  5,  5,  5,  5,  5,  6,  7,  8,  9, 10,
        /*G*/  7,  6,  5,  4,  4,  4,  5,  6,  6,  6,  6,  6,  6,  6,  7,  8,  9,
        /*T*/  8,  7,  6,  5,  5,  5,  5,  5,  6,  7,  7,  7,  7,  7,  7,  7,  8,
        /*A*/  9,  8,  7,  6,  6,  6,  6,  6,  6,  6,  7,  8,  8,  8,  8,  8,  8
        },
        std::vector
        {
        //     e,  A,  A,  C,  C,  G,  G,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/NON,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,
        /*A*/U  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,
        /*C*/U  ,U  ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL ,DL ,L  ,L  ,L  ,L  ,
        /*G*/U  ,U  ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL ,DL ,L  ,L  ,
        /*T*/U  ,U  ,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL ,DL ,
        /*A*/U  ,DU ,D  ,DUL,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,
        /*C*/U  ,U  ,U  ,D  ,DL ,DUL,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,
        /*G*/U  ,U  ,U  ,U  ,D  ,D  ,DL ,DUL,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,
        /*T*/U  ,U  ,U  ,U  ,DU ,DU ,D  ,D  ,DL ,DUL,DU ,DU ,DU ,DU ,D  ,D  ,DL ,
        /*A*/U  ,DU ,DU ,U  ,DU ,DU ,DU ,DU ,D  ,D  ,DL ,DUL,DU ,DU ,DU ,DU ,D
        }
    };
}();

static auto dna4_03 = []()
{
    return alignment_fixture
    {
        "TTTTTACGTATGTCCCCC"_dna4,
        "ACGTAAAACGT"_dna4,
        align_config_semi_seq2 | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
        -13,
        "-----ACGTA--------AAACGT",
        "TTTTTACGTATGTCCCCC------",
        alignment_coordinate{detail::column_index_type{0u}, detail::row_index_type{0u}},
        alignment_coordinate{detail::column_index_type{5u}, detail::row_index_type{18u}},
        std::vector
        {
        //     e,  A,  A,  C,  C,  G,  G,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
        /*A*/  1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
        /*C*/  2,  1,  1,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
        /*G*/  3,  2,  2,  2,  2,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13,
        /*T*/  4,  3,  3,  3,  3,  3,  3,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
        /*A*/  5,  4,  3,  4,  4,  4,  4,  4,  4,  4,  5,  6,  7,  8,  9, 10, 11,
        /*C*/  6,  5,  4,  3,  4,  5,  5,  5,  5,  5,  5,  5,  6,  7,  8,  9, 10,
        /*G*/  7,  6,  5,  4,  4,  4,  5,  6,  6,  6,  6,  6,  6,  6,  7,  8,  9,
        /*T*/  8,  7,  6,  5,  5,  5,  5,  5,  6,  7,  7,  7,  7,  7,  7,  7,  8,
        /*A*/  9,  8,  7,  6,  6,  6,  6,  6,  6,  6,  7,  8,  8,  8,  8,  8,  8
        },
        std::vector
        {
        //     e,  A,  A,  C,  C,  G,  G,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/NON,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,
        /*A*/U  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,
        /*C*/U  ,U  ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL ,DL ,L  ,L  ,L  ,L  ,
        /*G*/U  ,U  ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL ,DL ,L  ,L  ,
        /*T*/U  ,U  ,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL ,DL ,
        /*A*/U  ,DU ,D  ,DUL,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,
        /*C*/U  ,U  ,U  ,D  ,DL ,DUL,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,
        /*G*/U  ,U  ,U  ,U  ,D  ,D  ,DL ,DUL,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,
        /*T*/U  ,U  ,U  ,U  ,DU ,DU ,D  ,D  ,DL ,DUL,DU ,DU ,DU ,DU ,D  ,D  ,DL ,
        /*A*/U  ,DU ,DU ,U  ,DU ,DU ,DU ,DU ,D  ,D  ,DL ,DUL,DU ,DU ,DU ,DU ,D
        }
    };
}();

static auto dna4_04 = []()
{
    return alignment_fixture
    {
        "ACGTAAAACGT"_dna4,
        "TTTTTACGTATGTCCCCC"_dna4,
        align_config_semi_seq2 | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
        10,
        "TTTTTACGT---ATGTCCCCC",
        "-----ACGTAAAACGT-----",
        alignment_coordinate{detail::column_index_type{0u}, detail::row_index_type{0u}},
        alignment_coordinate{detail::column_index_type{13u}, detail::row_index_type{11u}},
        std::vector
        {
        //     e,  A,  A,  C,  C,  G,  G,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
        /*A*/  1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
        /*C*/  2,  1,  1,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
        /*G*/  3,  2,  2,  2,  2,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13,
        /*T*/  4,  3,  3,  3,  3,  3,  3,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
        /*A*/  5,  4,  3,  4,  4,  4,  4,  4,  4,  4,  5,  6,  7,  8,  9, 10, 11,
        /*C*/  6,  5,  4,  3,  4,  5,  5,  5,  5,  5,  5,  5,  6,  7,  8,  9, 10,
        /*G*/  7,  6,  5,  4,  4,  4,  5,  6,  6,  6,  6,  6,  6,  6,  7,  8,  9,
        /*T*/  8,  7,  6,  5,  5,  5,  5,  5,  6,  7,  7,  7,  7,  7,  7,  7,  8,
        /*A*/  9,  8,  7,  6,  6,  6,  6,  6,  6,  6,  7,  8,  8,  8,  8,  8,  8
        },
        std::vector
        {
        //     e,  A,  A,  C,  C,  G,  G,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/NON,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,L  ,
        /*A*/U  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,
        /*C*/U  ,U  ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL ,DL ,L  ,L  ,L  ,L  ,
        /*G*/U  ,U  ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL ,DL ,L  ,L  ,
        /*T*/U  ,U  ,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,DL ,DL ,
        /*A*/U  ,DU ,D  ,DUL,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,L  ,L  ,
        /*C*/U  ,U  ,U  ,D  ,DL ,DUL,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,L  ,L  ,
        /*G*/U  ,U  ,U  ,U  ,D  ,D  ,DL ,DUL,DU ,DU ,DU ,DU ,D  ,D  ,DL ,L  ,L  ,
        /*T*/U  ,U  ,U  ,U  ,DU ,DU ,D  ,D  ,DL ,DUL,DU ,DU ,DU ,DU ,D  ,D  ,DL ,
        /*A*/U  ,DU ,DU ,U  ,DU ,DU ,DU ,DU ,D  ,D  ,DL ,DUL,DU ,DU ,DU ,DU ,D
        }
    };
}();

} // namespace seqan3::test::alignment::fixture::global::affine::unbanded
