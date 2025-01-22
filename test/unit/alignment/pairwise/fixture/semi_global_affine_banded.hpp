// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "alignment_fixture.hpp"

using seqan3::operator""_dna4;

namespace seqan3::test::alignment::fixture::semi_global::affine::banded
{

inline constexpr auto config_band_m4_8 =
    seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-4}, seqan3::align_cfg::upper_diagonal{8}};

inline constexpr auto config_gap =
    seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10}, seqan3::align_cfg::extension_score{-1}};

inline constexpr auto config_semi_seq1 =
    seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                     seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                     seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                     seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}};

inline constexpr auto config_semi_seq2 =
    seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{false},
                                     seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                                     seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
                                     seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}};

// free gaps for the second leading and first trailing gaps.
inline constexpr auto config_free_gaps_sl_ft =
    seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{false},
                                     seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                                     seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                     seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}};

inline constexpr auto config_free_gaps_tlbr =
    seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                     seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                                     seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                     seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}};

inline constexpr auto config_scoring_m4_mm5 = seqan3::align_cfg::scoring_scheme{
    seqan3::nucleotide_scoring_scheme{seqan3::match_score{4}, seqan3::mismatch_score{-5}}};

static auto dna4_01_semi_first = []()
{
    return alignment_fixture{
        "TTTTTACGTATGTCCCCC"_dna4,
        "ACGTAAAACGT"_dna4,
        config_semi_seq1 | config_gap | config_band_m4_8 | config_scoring_m4_mm5,
        10,
        "ACGT---ATGT",
        "ACGTAAAACGT",
        /*.sequence1_begin_position = */ 5u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 13u,
        /*.sequence2_end_position = */ 11u,
        std::vector<std::optional<int32_t>>{
            //      e,  T,  T,  T,  T,  T,  A,  C,  G,  T,  A,  T,  G,  T,  C,  C,  C,  C,  C
            /*e*/ 0,   0,   0,   0,   0,   0,   0,   0,   0,  INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
            /*A*/ -11, -5,  -5,  -5,  -5,  -5,  4,   -5,  -5, -5,  INF, INF, INF, INF, INF, INF, INF, INF, INF,
            /*C*/ -12, -12, -10, -10, -10, -10, -7,  8,   -3, -4,  -5,  INF, INF, INF, INF, INF, INF, INF, INF,
            /*G*/ -13, -13, -13, -13, -13, -13, -8,  -3,  12, 1,   0,   -1,  INF, INF, INF, INF, INF, INF, INF,
            /*T*/ -14, -9,  -9,  -9,  -9,  -9,  -9,  -4,  1,  16,  5,   4,   3,   INF, INF, INF, INF, INF, INF,
            /*A*/ INF, -15, -14, -14, -14, -14, -5,  -5,  0,  5,   20,  9,   8,   7,   INF, INF, INF, INF, INF,
            /*A*/ INF, INF, -16, -16, -16, -16, -10, -6,  -1, 4,   9,   15,  4,   3,   2,   INF, INF, INF, INF,
            /*A*/ INF, INF, INF, -17, -17, -17, -12, -7,  -2, 3,   8,   4,   10,  -1,  -2,  -3,  INF, INF, INF,
            /*A*/ INF, INF, INF, INF, -18, -18, -13, -8,  -3, 2,   7,   3,   -1,  5,   -6,  -7,  -8,  INF, INF,
            /*C*/ INF, INF, INF, INF, INF, -19, -14, -9,  -4, 1,   6,   2,   -2,  -6,  9,   -2,  -3,  -4,  INF,
            /*G*/ INF, INF, INF, INF, INF, INF, -15, -10, -5, 0,   5,   1,   6,   -5,  -2,  4,   -7,  -8,  -9,
            /*T*/ INF, INF, INF, INF, INF, INF, INF, -11, -6, -1,  4,   9,   -2,  10,  -1,  -2,  -1,  -4,  -5},
        std::vector<std::optional<seqan3::detail::trace_directions>>{
            //      e,  T,  T,  T,  T,  T,  A,  C,  G,  T,  A,  T,  G,  T,  C,  C,  C,  C,  C
            /*e*/ N,   N,   N,   N,   N,   N,   N,   N,   N,   INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
            /*A*/ U,   DUL, DUL, DUL, DUL, DUL, DUL, DUL, DUl, D,   INF, INF, INF, INF, INF, INF, INF, INF, INF,
            /*C*/ u,   uL,  DuL, DuL, DuL, DuL, UL,  DuL, L,   l,   l,   INF, INF, INF, INF, INF, INF, INF, INF,
            /*G*/ u,   uL,  uL,  uL,  uL,  uL,  uL,  UL,  DuL, L,   l,   l,   INF, INF, INF, INF, INF, INF, INF,
            /*T*/ u,   DuL, DuL, DuL, DuL, DuL, uL,  uL,  UL,  DUL, L,   DUl, l,   INF, INF, INF, INF, INF, INF,
            /*A*/ INF, u,   DuL, DuL, DuL, DuL, DuL, uL,  uL,  UL,  DUL, L,   l,   l,   INF, INF, INF, INF, INF,
            /*A*/ INF, INF, u,   uL,  uL,  uL,  DuL, uL,  uL,  uL,  DUL, DUL, DUL, DUl, D,   INF, INF, INF, INF,
            /*A*/ INF, INF, INF, u,   uL,  uL,  DuL, uL,  uL,  uL,  DuL, DUL, Dul, DuL, DUl, D,   INF, INF, INF,
            /*A*/ INF, INF, INF, INF, u,   uL,  DuL, uL,  uL,  uL,  DuL, DuL, DUl, Dul, DuL, DUl, D,   INF, INF,
            /*C*/ INF, INF, INF, INF, INF, u,   uL,  DuL, uL,  uL,  uL,  DuL, Dul, DUl, Dul, DuL, DUl, D,   INF,
            /*G*/ INF, INF, INF, INF, INF, INF, u,   uL,  DuL, uL,  uL,  DuL, Dul, L,   Ul,  DUl, DUL, DUl, D,
            /*T*/ INF, INF, INF, INF, INF, INF, INF, u,   uL,  DuL, uL,  DuL, L,   Dul, L,   l,   Dul, l,   l}};
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
//         align_config_semi_seq1
//             | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
//                                                                                   seqan3::mismatch_score{-5}}},
//         -13,
//         "-----ACGTA--------AAACGT",
//         "TTTTTACGTATGTCCCCC------",
//         0u,
//         0u,
//         5u,
//         18u}
//     };
// }();

static auto dna4_03_semi_second = []()
{
    return alignment_fixture{
        "TTTTTACGTATGTCCCCC"_dna4,
        "ACGTAAAACGT"_dna4,
        config_semi_seq2 | config_gap | config_band_m4_8 | config_scoring_m4_mm5,
        -19,
        "TTTTTACGTATGTCCCCC",
        "GTAAAACGT---------",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 2u,
        /*.sequence1_end_position = */ 18u,
        /*.sequence2_end_position = */ 11u,
        std::vector<std::optional<int32_t>>{
            //      e,  T,  T,  T,  T,  T,  A,  C,  G,  T,  A,  T,  G,  T,  C,  C,  C,  C,  C
            /*e*/ 0,   -11, -12, -13, -14, -15, -16, -17, -18, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
            /*A*/ 0,   -5,  -12, -13, -14, -15, -11, -17, -18, -19, INF, INF, INF, INF, INF, INF, INF, INF, INF,
            /*C*/ 0,   -5,  -10, -13, -14, -15, -16, -7,  -18, -19, -20, INF, INF, INF, INF, INF, INF, INF, INF,
            /*G*/ 0,   -5,  -10, -13, -14, -15, -16, -17, -3,  -14, -15, -16, INF, INF, INF, INF, INF, INF, INF,
            /*T*/ 0,   4,   -1,  -6,  -9,  -10, -11, -12, -13, 1,   -10, -11, -12, INF, INF, INF, INF, INF, INF,
            /*A*/ INF, -5,  -1,  -6,  -11, -14, -6,  -16, -15, -10, 5,   -6,  -7,  -8,  INF, INF, INF, INF, INF,
            /*A*/ INF, INF, -10, -6,  -11, -16, -10, -11, -16, -11, -6,  0,   -11, -12, -13, INF, INF, INF, INF,
            /*A*/ INF, INF, INF, -15, -11, -16, -12, -15, -16, -12, -7,  -11, -5,  -16, -17, -18, INF, INF, INF,
            /*A*/ INF, INF, INF, INF, -20, -16, -12, -17, -18, -13, -8,  -12, -16, -10, -21, -22, -23, INF, INF,
            /*C*/ INF, INF, INF, INF, INF, -25, -20, -8,  -19, -14, -9,  -13, -17, -21, -6,  -17, -18, -19, INF,
            /*G*/ INF, INF, INF, INF, INF, INF, -21, -19, -4,  -15, -10, -14, -9,  -19, -17, -11, -22, -23, -24,
            /*T*/ INF, INF, INF, INF, INF, INF, INF, -20, -15, 0,   -11, -6,  -13, -5,  -15, -16, -16, -18, -19},
        std::vector<std::optional<seqan3::detail::trace_directions>>{
            //      e,  T,  T,  T,  T,  T,  A,  C,  G,  T,  A,  T,  G,  T,  C,  C,  C,  C,  C
            /*e*/ N,   L,   l,   l,   l,   l,   l,   l,   l,   INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
            /*A*/ N,   DUL, l,   l,   l,   l,   DUl, l,   l,   l,   INF, INF, INF, INF, INF, INF, INF, INF, INF,
            /*C*/ N,   DUL, DUl, l,   l,   l,   l,   DUl, l,   l,   l,   INF, INF, INF, INF, INF, INF, INF, INF,
            /*G*/ N,   DUL, DUl, l,   l,   l,   l,   l,   DUl, L,   l,   l,   INF, INF, INF, INF, INF, INF, INF,
            /*T*/ N,   DUL, DUL, DUl, DUl, DUl, l,   l,   l,   DUl, L,   DUl, l,   INF, INF, INF, INF, INF, INF,
            /*A*/ INF, DU,  DUL, DUL, DUl, DUl, DUl, Dul, ul,  Ul,  DUl, L,   l,   l,   INF, INF, INF, INF, INF,
            /*A*/ INF, INF, DU,  DUL, DuL, Dul, DUl, Dul, ul,  ul,  DUL, DUL, DUL, DUl, D,   INF, INF, INF, INF,
            /*A*/ INF, INF, INF, DU,  DuL, DuL, Dul, DuL, Dul, ul,  DuL, DUL, Dul, DuL, DUl, D,   INF, INF, INF,
            /*A*/ INF, INF, INF, INF, DU,  DuL, DuL, DuL, ul,  ul,  DuL, DuL, DUl, Dul, DuL, DUl, D,   INF, INF,
            /*C*/ INF, INF, INF, INF, INF, Du,  uL,  DuL, uL,  ul,  ul,  DuL, Dul, DUl, Dul, DuL, DUl, D,   INF,
            /*G*/ INF, INF, INF, INF, INF, INF, u,   UL,  DuL, uL,  ul,  Dul, Dul, l,   Ul,  DUl, DUl, DUl, D,
            /*T*/ INF, INF, INF, INF, INF, INF, INF, u,   UL,  DuL, uL,  Dul, l,   Dul, l,   l,   Dul, l,   l}};
}();

static auto dna4_04_semi_second = []()
{
    return alignment_fixture{
        "ACGTAAAACGT"_dna4,
        "TTTTTACGTATGTCCCCC"_dna4,
        config_semi_seq2 | config_gap | config_band_m4_8 | config_scoring_m4_mm5,
        -5,
        "ACGTAAAACGT",
        "------TACGT",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 4u,
        /*.sequence1_end_position = */ 11u,
        /*.sequence2_end_position = */ 9u,
        std::vector<std::optional<int32_t>>{//      e,  A,  C,  G,  T,  A,  A,  A,  A,  C,  G,  T
                                            /*e*/ 0,   -11, -12, -13, -14, -15, -16, -17, -18, INF, INF, INF,
                                            /*T*/ 0,   -5,  -12, -13, -9,  -15, -16, -17, -18, -19, INF, INF,
                                            /*T*/ 0,   -5,  -10, -13, -9,  -14, -16, -17, -18, -19, -20, INF,
                                            /*T*/ 0,   -5,  -10, -13, -9,  -14, -16, -17, -18, -19, -20, -16,
                                            /*T*/ 0,   -5,  -10, -13, -9,  -14, -16, -17, -18, -19, -20, -16,
                                            /*T*/ INF, -5,  -10, -15, -9,  -14, -19, -21, -22, -23, -24, -16,
                                            /*A*/ INF, INF, -10, -15, -20, -5,  -10, -15, -17, -19, -20, -21,
                                            /*C*/ INF, INF, INF, -15, -20, -16, -10, -15, -20, -13, -24, -25,
                                            /*G*/ INF, INF, INF, INF, -20, -17, -21, -15, -20, -24, -9,  -20,
                                            /*T*/ INF, INF, INF, INF, INF, -18, -22, -26, -20, -25, -20, -5,
                                            /*A*/ INF, INF, INF, INF, INF, INF, -14, -18, -22, -25, -21, -16,
                                            /*T*/ INF, INF, INF, INF, INF, INF, INF, -19, -23, -27, -22, -17,
                                            /*G*/ INF, INF, INF, INF, INF, INF, INF, INF, -24, -28, -23, -18,
                                            /*T*/ INF, INF, INF, INF, INF, INF, INF, INF, INF, -29, -24, -19,
                                            /*C*/ INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, -25, -20,
                                            /*C*/ INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, -21,
                                            /*C*/ INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
                                            /*C*/ INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
                                            /*C*/ INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF},
        std::vector<std::optional<seqan3::detail::trace_directions>>{
            //      e,  A,  C,  G,  T,  A,  A,  A,  A,  C,  G,  T
            /*e*/ N,   L,   l,   l,   l,   l,   l,   l,   l,   INF, INF, INF,
            /*T*/ N,   DUL, l,   l,   DUl, l,   l,   l,   l,   l,   INF, INF,
            /*T*/ N,   DUL, DUl, l,   DUl, DUl, l,   l,   l,   l,   l,   INF,
            /*T*/ N,   DUL, DUl, l,   DUl, DUl, l,   l,   l,   l,   l,   D,
            /*T*/ N,   DUL, DUl, l,   DUl, DUl, l,   l,   l,   l,   l,   DUl,
            /*T*/ INF, DU,  DUL, DUl, DUl, DUl, DUl, DUl, DUl, DUl, DUl, DUl,
            /*A*/ INF, INF, DU,  DuL, DUl, DUl, DuL, Dul, Dul, l,   l,   l,
            /*C*/ INF, INF, INF, Du,  DuL, Ul,  DUL, DUL, DUl, DUl, DUl, Dul,
            /*G*/ INF, INF, INF, INF, Du,  uL,  DUL, DUl, DuL, Ul,  Dul, L,
            /*T*/ INF, INF, INF, INF, INF, u,   DuL, DUl, Dul, DuL, Ul,  DuL,
            /*A*/ INF, INF, INF, INF, INF, INF, Du,  DuL, Dul, Dul, ul,  Ul,
            /*T*/ INF, INF, INF, INF, INF, INF, INF, Du,  DuL, Dul, ul,  Dul,
            /*G*/ INF, INF, INF, INF, INF, INF, INF, INF, Du,  DuL, Dul, uL,
            /*T*/ INF, INF, INF, INF, INF, INF, INF, INF, INF, Du,  uL,  DuL,
            /*C*/ INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, u,   uL,
            /*C*/ INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, u,
            /*C*/ INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
            /*C*/ INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
            /*C*/ INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF}};
}();

static auto dna4_free_lb_with_band_tl2br_no_matches = []()
{
    return alignment_fixture{
        "AAAAAAAAAAA"_dna4,
        "TTTTTTTTTTT"_dna4,
        config_free_gaps_sl_ft | config_gap | config_band_m4_8 | config_scoring_m4_mm5,
        -34,
        "AAAAAAA-------",
        "-------TTTTTTT",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 4u,
        /*.sequence1_end_position = */ 7u,
        /*.sequence2_end_position = */ 11u,
        std::vector<std::optional<int32_t>>{//      e,  A,  A,  A,  A,  A,  A,  A,  A,  A,  A,  A,
                                            /*e*/ 0,   -11, -12, -13, -14, -15, -16, -17, -18, INF, INF, INF,
                                            /*T*/ 0,   -5,  -12, -13, -14, -15, -16, -17, -18, -19, INF, INF,
                                            /*T*/ 0,   -5,  -10, -13, -14, -15, -16, -17, -18, -19, -20, INF,
                                            /*T*/ 0,   -5,  -10, -13, -14, -15, -16, -17, -18, -19, -20, -21,
                                            /*T*/ 0,   -5,  -10, -13, -14, -15, -16, -17, -18, -19, -20, -21,
                                            /*T*/ INF, -5,  -10, -15, -18, -19, -20, -21, -22, -23, -24, -25,
                                            /*T*/ INF, INF, -10, -15, -20, -23, -24, -25, -26, -27, -28, -29,
                                            /*T*/ INF, INF, INF, -15, -20, -25, -28, -29, -30, -31, -32, -33,
                                            /*T*/ INF, INF, INF, INF, -20, -25, -30, -31, -32, -33, -34, -35,
                                            /*T*/ INF, INF, INF, INF, INF, -25, -30, -32, -33, -34, -35, -36,
                                            /*T*/ INF, INF, INF, INF, INF, INF, -30, -33, -34, -35, -36, -37,
                                            /*T*/ INF, INF, INF, INF, INF, INF, INF, -34, -35, -36, -37, -38},
        std::vector<std::optional<seqan3::detail::trace_directions>>{
            //      e,  A,  A,  A,  A,  A,  A,  A,  A,  A,  A,  A,
            /*e*/ N,   L,   l,   l,   l,   l,   l,   l,   l,   INF, INF, INF,
            /*T*/ N,   DUL, l,   l,   l,   l,   l,   l,   l,   l,   INF, INF,
            /*T*/ N,   DUL, DUl, l,   l,   l,   l,   l,   l,   l,   l,   INF,
            /*T*/ N,   DUL, DUl, l,   l,   l,   l,   l,   l,   l,   l,   l,
            /*T*/ N,   DUL, DUl, l,   l,   l,   l,   l,   l,   l,   l,   l,
            /*T*/ INF, DU,  DUL, DUl, DUl, DUl, DUl, DUl, DUl, DUl, DUl, DUl,
            /*T*/ INF, INF, DU,  DuL, Dul, Dul, Dul, Dul, Dul, Dul, Dul, Dul,
            /*T*/ INF, INF, INF, Du,  DuL, Dul, Dul, Dul, Dul, Dul, Dul, Dul,
            /*T*/ INF, INF, INF, INF, Du,  DuL, Dul, ul,  ul,  ul,  ul,  ul,
            /*T*/ INF, INF, INF, INF, INF, Du,  DuL, ul,  ul,  ul,  ul,  ul,
            /*T*/ INF, INF, INF, INF, INF, INF, Du,  uL,  ul,  ul,  ul,  ul,
            /*T*/ INF, INF, INF, INF, INF, INF, INF, u,   uL,  ul,  ul,  ul}};
}();

// band starts in top left with lower diagonal exceeding size of sequence 2 and ends in the bottom.
static auto dna4_free_tlbr_with_band_tl2b = []()
{
    return alignment_fixture{
        "AGATTTACTACGCAT"_dna4,
        "GTAGCAT"_dna4,
        config_free_gaps_tlbr | config_gap | config_scoring_m4_mm5
            | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-10},
                                                 seqan3::align_cfg::upper_diagonal{4}},
        5,
        "AG-AT",
        "AGCAT",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 2u,
        /*.sequence1_end_position = */ 4u,
        /*.sequence2_end_position = */ 7u,
        std::vector<std::optional<int32_t>>{
            //    e  ,A  ,G  ,A  ,T  ,T  ,T  ,A  ,C  ,T  ,A  ,C  ,G  ,C  ,A  ,T  ,
            /*e*/ 0, 0,  0,  0,  0,  INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
            /*G*/ 0, -5, 4,  -5, -5, -5,  INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
            /*T*/ 0, -5, -7, -1, -1, -1,  -1,  INF, INF, INF, INF, INF, INF, INF, INF, INF,
            /*A*/ 0, 4,  -7, -3, -6, -6,  -6,  3,   INF, INF, INF, INF, INF, INF, INF, INF,
            /*G*/ 0, -5, 8,  -3, -4, -5,  -6,  -7,  -2,  INF, INF, INF, INF, INF, INF, INF,
            /*C*/ 0, -5, -3, 3,  -8, -9,  -10, -9,  -3,  -7,  INF, INF, INF, INF, INF, INF,
            /*A*/ 0, 4,  -4, 1,  -2, -10, -11, -6,  -13, -8,  -3,  INF, INF, INF, INF, INF,
            /*T*/ 0, -5, -1, -9, 5,  2,   -6,  -8,  -9,  -9,  -11, -8,  INF, INF, INF, INF},
        std::vector<std::optional<seqan3::detail::trace_directions>>{
            //    e  ,A  ,G  ,A  ,T  ,T  ,T  ,A  ,C  ,T  ,A  ,C  ,G  ,C  ,A  ,T  ,
            /*e*/ N, N,   N,   N,   N,   INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
            /*G*/ N, DUL, DUl, DUL, DUl, D,   INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
            /*T*/ N, DuL, Ul,  Dul, DuL, DUL, D,   INF, INF, INF, INF, INF, INF, INF, INF, INF,
            /*A*/ N, DuL, L,   DUl, DUl, DUl, DUl, D,   INF, INF, INF, INF, INF, INF, INF, INF,
            /*G*/ N, DUL, Dul, L,   l,   l,   l,   l,   D,   INF, INF, INF, INF, INF, INF, INF,
            /*C*/ N, DuL, Ul,  Dul, DuL, Dul, Dul, ul,  DUl, D,   INF, INF, INF, INF, INF, INF,
            /*A*/ N, DuL, uL,  DUl, Dul, l,   l,   Dul, l,   DUl, D,   INF, INF, INF, INF, INF,
            /*T*/ N, DUL, Dul, DuL, DUl, DuL, Dul, l,   l,   Dul, l,   D,   INF, INF, INF, INF}};
}();

// band starts in top left with upper diagonal exceeding size of sequence 1 and ends in the right side.
static auto dna4_free_tlbr_with_band_tl2r = []()
{
    return alignment_fixture{"GTAGCAT"_dna4,
                             "AGATTTACTACGCAT"_dna4,
                             config_free_gaps_tlbr | config_gap | config_scoring_m4_mm5
                                 | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-3},
                                                                      seqan3::align_cfg::upper_diagonal{100}},
                             5,
                             "AGCAT",
                             "AG-AT",
                             /*.sequence1_begin_position = */ 2u,
                             /*.sequence2_begin_position = */ 0u,
                             /*.sequence1_end_position = */ 7u,
                             /*.sequence2_end_position = */ 4u,
                             std::vector<std::optional<int32_t>>{//    e  ,G  ,T  ,A  ,G  ,C  ,A  ,T  ,
                                                                 /*e*/ 0,   0,   0,   0,   0,   0,   0,   0,
                                                                 /*A*/ 0,   -5,  -5,  4,   -5,  -5,  4,   -5,
                                                                 /*G*/ 0,   4,   -7,  -7,  8,   -3,  -4,  -1,
                                                                 /*A*/ 0,   -5,  -1,  -3,  -3,  3,   1,   -9,
                                                                 /*T*/ INF, -5,  -1,  -6,  -4,  -8,  -2,  5,
                                                                 /*T*/ INF, INF, -1,  -6,  -5,  -9,  -10, 2,
                                                                 /*T*/ INF, INF, INF, -6,  -6,  -10, -11, -6,
                                                                 /*A*/ INF, INF, INF, INF, -7,  -11, -6,  -8,
                                                                 /*C*/ INF, INF, INF, INF, INF, -3,  -13, -9,
                                                                 /*T*/ INF, INF, INF, INF, INF, INF, -8,  -9,
                                                                 /*A*/ INF, INF, INF, INF, INF, INF, INF, -11,
                                                                 /*C*/ INF, INF, INF, INF, INF, INF, INF, INF,
                                                                 /*G*/ INF, INF, INF, INF, INF, INF, INF, INF,
                                                                 /*C*/ INF, INF, INF, INF, INF, INF, INF, INF,
                                                                 /*A*/ INF, INF, INF, INF, INF, INF, INF, INF,
                                                                 /*T*/ INF, INF, INF, INF, INF, INF, INF, INF},
                             std::vector<std::optional<seqan3::detail::trace_directions>>{
                                 //    e  ,G  ,T  ,A  ,G  ,C  ,A  ,T  ,
                                 /*e*/ N,   N,   N,   N,   N,   N,   N,   N,
                                 /*A*/ N,   DUL, DUl, DUl, DUL, DUl, DUl, DUL,
                                 /*G*/ N,   DuL, L,   Ul,  Dul, L,   l,   Dul,
                                 /*A*/ N,   DUL, Dul, DuL, Ul,  Dul, DuL, DUl,
                                 /*T*/ INF, Du,  DUL, DuL, ul,  DUl, Dul, DuL,
                                 /*T*/ INF, INF, DU,  DuL, ul,  Dul, ul,  DUl,
                                 /*T*/ INF, INF, INF, Du,  uL,  DuL, ul,  Dul,
                                 /*A*/ INF, INF, INF, INF, u,   DuL, Dul, uL,
                                 /*C*/ INF, INF, INF, INF, INF, Du,  uL,  ul,
                                 /*T*/ INF, INF, INF, INF, INF, INF, Du,  DuL,
                                 /*A*/ INF, INF, INF, INF, INF, INF, INF, u,
                                 /*C*/ INF, INF, INF, INF, INF, INF, INF, INF,
                                 /*G*/ INF, INF, INF, INF, INF, INF, INF, INF,
                                 /*C*/ INF, INF, INF, INF, INF, INF, INF, INF,
                                 /*A*/ INF, INF, INF, INF, INF, INF, INF, INF,
                                 /*T*/ INF, INF, INF, INF, INF, INF, INF, INF}};
}();

} // namespace seqan3::test::alignment::fixture::semi_global::affine::banded
