// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

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

namespace seqan3::test::alignment::fixture::semi_global::affine::unbanded
{
inline constexpr auto align_config =
    seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10}, seqan3::align_cfg::extension_score{-1}};

inline constexpr auto align_config_semi_seq1 =
    align_config
    | seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                       seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                       seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                       seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}};
inline constexpr auto align_config_semi_seq2 =
    align_config
    | seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{false},
                                       seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                                       seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
                                       seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}};

static auto dna4_01_semi_first = []()
{
    using seqan3::detail::column_index_type;
    using seqan3::detail::row_index_type;

    return alignment_fixture{
        "TTTTTACGTATGTCCCCC"_dna4,
        "ACGTAAAACGT"_dna4,
        align_config_semi_seq1
            | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                  seqan3::mismatch_score{-5}}},
        10,
        "ACGT---ATGT",
        "ACGTAAAACGT",
        /*.sequence1_begin_position = */ 5u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 13u,
        /*.sequence2_end_position = */ 11u,
        std::vector{//      e,  T,  T,  T,  T,  T,  A,  C,  G,  T,  A,  T,  G,  T,  C,  C,  C,  C,  C
                    /*e*/ 0,   0,   0,   0,   0,   0,   0,   0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                    /*A*/ -11, -5,  -5,  -5,  -5,  -5,  4,   -5,  -5, -5, 4,  -5, -5, -5, -5, -5, -5, -5, -5,
                    /*C*/ -12, -12, -10, -10, -10, -10, -7,  8,   -3, -4, -5, -1, -7, -8, -1, -1, -1, -1, -1,
                    /*G*/ -13, -13, -13, -13, -13, -13, -8,  -3,  12, 1,  0,  -1, 3,  -3, -4, -5, -6, -6, -6,
                    /*T*/ -14, -9,  -9,  -9,  -9,  -9,  -9,  -4,  1,  16, 5,  4,  3,  7,  1,  0,  -1, -2, -3,
                    /*A*/ -15, -15, -14, -14, -14, -14, -5,  -5,  0,  5,  20, 9,  8,  7,  6,  5,  4,  3,  2,
                    /*A*/ -16, -16, -16, -16, -16, -16, -10, -6,  -1, 4,  9,  15, 4,  3,  2,  1,  0,  -1, -2,
                    /*A*/ -17, -17, -17, -17, -17, -17, -12, -7,  -2, 3,  8,  4,  10, -1, -2, -3, -4, -5, -6,
                    /*A*/ -18, -18, -18, -18, -18, -18, -13, -8,  -3, 2,  7,  3,  -1, 5,  -6, -7, -8, -9, -10,
                    /*C*/ -19, -19, -19, -19, -19, -19, -14, -9,  -4, 1,  6,  2,  -2, -6, 9,  -2, -3, -4, -5,
                    /*G*/ -20, -20, -20, -20, -20, -20, -15, -10, -5, 0,  5,  1,  6,  -5, -2, 4,  -7, -8, -9,
                    /*T*/ -21, -16, -16, -16, -16, -16, -16, -11, -6, -1, 4,  9,  -2, 10, -1, -2, -1, -4, -5},
        std::vector{
            //      e,  T,  T,  T,  T,  T,  A,  C,  G,  T,  A,  T,  G,  T,  C,  C,  C,  C,  C
            /*e*/ N, N,   N,   N,   N,   N,   N,   N,   N,   N,   N,   N,   N,   N,   N,   N,   N,   N,   N,
            /*A*/ U, DUL, DUL, DUL, DUL, DUL, DUL, DUL, DUl, DUl, DUl, DUL, DUl, DUl, DUl, DUl, DUl, DUl, DUl,
            /*C*/ u, uL,  DuL, DuL, DuL, DuL, UL,  DuL, L,   l,   l,   Dul, l,   l,   Dul, Dul, Dul, Dul, DuL,
            /*G*/ u, uL,  uL,  uL,  uL,  uL,  uL,  UL,  DuL, L,   l,   l,   Dul, l,   l,   l,   DUl, DUl, DUl,
            /*T*/ u, DuL, DuL, DuL, DuL, DuL, uL,  uL,  UL,  DUL, L,   DUl, l,   Dul, l,   l,   l,   l,   l,
            /*A*/ u, uL,  DuL, DuL, DuL, DuL, DuL, uL,  uL,  UL,  DUL, L,   l,   l,   l,   l,   l,   l,   l,
            /*A*/ u, uL,  uL,  uL,  uL,  uL,  DuL, uL,  uL,  uL,  DUL, DUL, DUL, DUl, DUl, DUl, DUl, DUl, DUl,
            /*A*/ u, uL,  uL,  uL,  uL,  uL,  DuL, uL,  uL,  uL,  DuL, DUL, Dul, DuL, Dul, Dul, Dul, Dul, Dul,
            /*A*/ u, uL,  uL,  uL,  uL,  uL,  DuL, uL,  uL,  uL,  DuL, DuL, DUl, Dul, DuL, Dul, Dul, Dul, Dul,
            /*C*/ u, uL,  uL,  uL,  uL,  uL,  uL,  DuL, uL,  uL,  uL,  DuL, Dul, DUl, Dul, DuL, Dul, Dul, Dul,
            /*G*/ u, uL,  uL,  uL,  uL,  uL,  uL,  uL,  DuL, uL,  uL,  DuL, Dul, L,   Ul,  Dul, DuL, Dul, Dul,
            /*T*/ u, DuL, DuL, DuL, DuL, DuL, uL,  uL,  uL,  DuL, uL,  DuL, L,   Dul, L,   l,   Dul, l,   l}};
}();

static auto dna4_02_semi_first = []()
{
    using seqan3::detail::column_index_type;
    using seqan3::detail::row_index_type;

    return alignment_fixture{
        "ACGTAAAACGT"_dna4,
        "TTTTTACGTATGTCCCCC"_dna4,
        align_config_semi_seq1
            | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                  seqan3::mismatch_score{-5}}},
        -13,
        "-----ACGTA--------",
        "TTTTTACGTATGTCCCCC",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 5u,
        /*.sequence2_end_position = */ 18u,
        std::vector{//      e,  A,  C,  G,  T,  A,  A,  A,  A,  C,  G,  T
                    /*e*/ 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
                    /*T*/ -11, -5,  -5,  -5,  4,   -5,  -5,  -5,  -5,  -5,  -5,  4,
                    /*T*/ -12, -12, -10, -10, -1,  -1,  -10, -10, -10, -10, -10, -1,
                    /*T*/ -13, -13, -13, -13, -6,  -6,  -6,  -13, -13, -13, -13, -6,
                    /*T*/ -14, -14, -14, -14, -9,  -11, -11, -11, -14, -14, -14, -9,
                    /*T*/ -15, -15, -15, -15, -10, -14, -15, -15, -15, -15, -15, -10,
                    /*A*/ -16, -11, -16, -16, -11, -6,  -10, -11, -11, -16, -16, -11,
                    /*C*/ -17, -17, -7,  -17, -12, -16, -11, -15, -16, -7,  -17, -12,
                    /*G*/ -18, -18, -18, -3,  -13, -15, -16, -16, -18, -18, -3,  -13,
                    /*T*/ -19, -19, -19, -14, 1,   -10, -11, -12, -13, -14, -14, 1,
                    /*A*/ -20, -15, -20, -15, -10, 5,   -6,  -7,  -8,  -9,  -10, -10,
                    /*T*/ -21, -21, -20, -16, -11, -6,  0,   -11, -12, -13, -14, -6,
                    /*G*/ -22, -22, -22, -16, -12, -7,  -11, -5,  -16, -17, -9,  -12,
                    /*T*/ -23, -23, -23, -18, -12, -8,  -12, -16, -10, -21, -18, -5,
                    /*C*/ -24, -24, -19, -19, -14, -9,  -13, -17, -21, -6,  -17, -14,
                    /*C*/ -25, -25, -20, -20, -15, -10, -14, -18, -22, -17, -11, -15,
                    /*C*/ -26, -26, -21, -21, -16, -11, -15, -19, -23, -18, -21, -16,
                    /*C*/ -27, -27, -22, -22, -17, -12, -16, -20, -24, -19, -22, -17,
                    /*C*/ -28, -28, -23, -23, -18, -13, -17, -21, -25, -20, -23, -18},
        std::vector{//      e,  A,  C,  G,  T,  A,  A,  A,  A,  C,  G,  T
                    /*e*/ N, N,   N,   N,   N,   N,   N,   N,   N,   N,   N,   N,
                    /*T*/ U, DUL, DUL, DUL, DUL, DUL, DUl, DUl, DUl, DUl, DUl, DUl,
                    /*T*/ u, uL,  DuL, DuL, DUL, DuL, DuL, Dul, Dul, Dul, Dul, DUl,
                    /*T*/ u, uL,  uL,  uL,  DuL, DUL, DuL, uL,  ul,  ul,  ul,  Dul,
                    /*T*/ u, uL,  uL,  uL,  DuL, DuL, Dul, Dul, uL,  ul,  ul,  Dul,
                    /*T*/ u, uL,  uL,  uL,  DuL, DuL, ul,  ul,  ul,  ul,  ul,  DuL,
                    /*A*/ u, DuL, uL,  ul,  ul,  DuL, DuL, Dul, Dul, ul,  ul,  ul,
                    /*C*/ u, uL,  DuL, uL,  ul,  Dul, Dul, Dul, Dul, Dul, uL,  ul,
                    /*G*/ u, uL,  uL,  DuL, uL,  l,   l,   Dul, ul,  ul,  Dul, uL,
                    /*T*/ u, uL,  uL,  UL,  DuL, L,   l,   l,   l,   l,   Ul,  Dul,
                    /*A*/ u, DuL, uL,  ul,  UL,  DuL, DuL, Dul, Dul, l,   l,   Ul,
                    /*T*/ u, uL,  DuL, uL,  DuL, UL,  DUL, DUL, DUl, DUl, Dul, Dul,
                    /*G*/ u, uL,  uL,  DuL, uL,  uL,  DUL, Dul, DuL, Dul, Dul, ul,
                    /*T*/ u, uL,  uL,  uL,  DuL, uL,  DuL, DUl, Dul, DuL, ul,  Dul,
                    /*C*/ u, uL,  DuL, uL,  uL,  uL,  DuL, Dul, DUl, Dul, L,   ul,
                    /*C*/ u, uL,  DuL, uL,  uL,  uL,  DuL, Dul, Dul, DUl, Dul, uL,
                    /*C*/ u, uL,  DuL, uL,  uL,  uL,  DuL, Dul, Dul, Dul, ul,  Dul,
                    /*C*/ u, uL,  DuL, uL,  uL,  uL,  DuL, Dul, Dul, Dul, ul,  ul,
                    /*C*/ u, uL,  DuL, uL,  uL,  uL,  DuL, Dul, Dul, Dul, ul,  ul}};
}();

static auto dna4_03_semi_second = []()
{
    using seqan3::detail::column_index_type;
    using seqan3::detail::row_index_type;

    return alignment_fixture{
        "TTTTTACGTATGTCCCCC"_dna4,
        "ACGTAAAACGT"_dna4,
        align_config_semi_seq2
            | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                  seqan3::mismatch_score{-5}}},
        -13,
        "TTTTTACGTATGTCCCCC",
        "-----ACGTA--------",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 18u,
        /*.sequence2_end_position = */ 5u,
        std::vector{
            //      e,  T,  T,  T,  T,  T,  A,  C,  G,  T,  A,  T,  G,  T,  C,  C,  C,  C,  C
            /*e*/ 0, -11, -12, -13, -14, -15, -16, -17, -18, -19, -20, -21, -22, -23, -24, -25, -26, -27, -28,
            /*A*/ 0, -5,  -12, -13, -14, -15, -11, -17, -18, -19, -15, -21, -22, -23, -24, -25, -26, -27, -28,
            /*C*/ 0, -5,  -10, -13, -14, -15, -16, -7,  -18, -19, -20, -20, -22, -23, -19, -20, -21, -22, -23,
            /*G*/ 0, -5,  -10, -13, -14, -15, -16, -17, -3,  -14, -15, -16, -16, -18, -19, -20, -21, -22, -23,
            /*T*/ 0, 4,   -1,  -6,  -9,  -10, -11, -12, -13, 1,   -10, -11, -12, -12, -14, -15, -16, -17, -18,
            /*A*/ 0, -5,  -1,  -6,  -11, -14, -6,  -16, -15, -10, 5,   -6,  -7,  -8,  -9,  -10, -11, -12, -13,
            /*A*/ 0, -5,  -10, -6,  -11, -15, -10, -11, -16, -11, -6,  0,   -11, -12, -13, -14, -15, -16, -17,
            /*A*/ 0, -5,  -10, -13, -11, -15, -11, -15, -16, -12, -7,  -11, -5,  -16, -17, -18, -19, -20, -21,
            /*A*/ 0, -5,  -10, -13, -14, -15, -11, -16, -18, -13, -8,  -12, -16, -10, -21, -22, -23, -24, -25,
            /*C*/ 0, -5,  -10, -13, -14, -15, -16, -7,  -18, -14, -9,  -13, -17, -21, -6,  -17, -18, -19, -20,
            /*G*/ 0, -5,  -10, -13, -14, -15, -16, -17, -3,  -14, -10, -14, -9,  -18, -17, -11, -21, -22, -23,
            /*T*/ 0, 4,   -1,  -6,  -9,  -10, -11, -12, -13, 1,   -10, -6,  -12, -5,  -14, -15, -16, -17, -18},
        std::vector{
            //      e,  T,  T,  T,  T,  T,  A,  C,  G,  T,  A,  T,  G,  T,  C,  C,  C,  C,  C
            /*e*/ N, L,   l,   l,   l,   l,   l,   l,   l,   l,   l,   l,   l,   l,   l,   l,   l,   l,   l,
            /*A*/ N, DUL, l,   l,   l,   l,   DUl, l,   l,   l,   DUl, l,   l,   l,   l,   l,   l,   l,   l,
            /*C*/ N, DUL, DUl, l,   l,   l,   l,   DUl, l,   l,   l,   DUl, l,   l,   DUl, DUl, DUl, DUl, DUl,
            /*G*/ N, DUL, DUl, l,   l,   l,   l,   l,   DUl, L,   l,   l,   DUl, l,   l,   l,   l,   l,   l,
            /*T*/ N, DUL, DUL, DUl, DUl, DUl, l,   l,   l,   DUl, L,   DUl, l,   DUl, l,   l,   l,   l,   l,
            /*A*/ N, DUL, DUl, DUL, DUl, DUl, DUl, Dul, ul,  Ul,  DUl, L,   l,   l,   l,   l,   l,   l,   l,
            /*A*/ N, DuL, DUl, DUl, Dul, l,   DUl, Dul, ul,  ul,  DUl, DUL, DUL, DUl, DUl, DUl, DUl, DUl, DUl,
            /*A*/ N, DuL, Dul, l,   Dul, l,   Dul, Dul, Dul, ul,  Dul, DUL, Dul, DuL, Dul, Dul, Dul, Dul, Dul,
            /*A*/ N, DuL, Dul, l,   l,   l,   Dul, Dul, ul,  ul,  Dul, DuL, DUl, Dul, DuL, Dul, Dul, Dul, Dul,
            /*C*/ N, DuL, Dul, l,   l,   l,   l,   Dul, l,   ul,  ul,  DuL, Dul, DUl, Dul, DuL, Dul, Dul, Dul,
            /*G*/ N, DuL, Dul, l,   l,   l,   l,   l,   Dul, L,   ul,  Dul, Dul, l,   Ul,  Dul, l,   l,   l,
            /*T*/ N, DuL, DuL, Dul, Dul, DUl, l,   l,   l,   Dul, L,   Dul, l,   Dul, l,   l,   Dul, l,   l}};
}();

static auto dna4_04_semi_second = []()
{
    using seqan3::detail::column_index_type;
    using seqan3::detail::row_index_type;

    return alignment_fixture{
        "ACGTAAAACGT"_dna4,
        "TTTTTACGTATGTCCCCC"_dna4,
        align_config_semi_seq2
            | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                  seqan3::mismatch_score{-5}}},
        10,
        "ACGTAAAACGT",
        "ACGT---ATGT",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 5u,
        /*.sequence1_end_position = */ 11u,
        /*.sequence2_end_position = */ 13u,
        std::vector{
            //      e,  A,  C,  G,  T,  A,  A,  A,  A,  C,  G,  T
            /*e*/ 0, -11, -12, -13, -14, -15, -16, -17, -18, -19, -20, -21,
            /*T*/ 0, -5,  -12, -13, -9,  -15, -16, -17, -18, -19, -20, -16,
            /*T*/ 0, -5,  -10, -13, -9,  -14, -16, -17, -18, -19, -20, -16,
            /*T*/ 0, -5,  -10, -13, -9,  -14, -16, -17, -18, -19, -20, -16,
            /*T*/ 0, -5,  -10, -13, -9,  -14, -16, -17, -18, -19, -20, -16,
            /*T*/ 0, -5,  -10, -13, -9,  -14, -16, -17, -18, -19, -20, -16,
            /*A*/ 0, 4,   -7,  -8,  -9,  -5,  -10, -12, -13, -14, -15, -16,
            /*C*/ 0, -5,  8,   -3,  -4,  -5,  -6,  -7,  -8,  -9,  -10, -11,
            /*G*/ 0, -5,  -3,  12,  1,   0,   -1,  -2,  -3,  -4,  -5,  -6,
            /*T*/ 0, -5,  -4,  1,   16,  5,   4,   3,   2,   1,   0,   -1,
            /*A*/ 0, 4,   -5,  0,   5,   20,  9,   8,   7,   6,   5,   4,
            /*T*/ 0, -5,  -1,  -1,  4,   9,   15,  4,   3,   2,   1,   9,
            /*G*/ 0, -5,  -7,  3,   3,   8,   4,   10,  -1,  -2,  6,   -2,
            /*T*/ 0, -5,  -8,  -3,  7,   7,   3,   -1,  5,   -6,  -5,  10,
            /*C*/ 0, -5,  -1,  -4,  1,   6,   2,   -2,  -6,  9,   -2,  -1,
            /*C*/ 0, -5,  -1,  -5,  0,   5,   1,   -3,  -7,  -2,  4,   -2,
            /*C*/ 0, -5,  -1,  -6,  -1,  4,   0,   -4,  -8,  -3,  -7,  -1,
            /*C*/ 0, -5,  -1,  -6,  -2,  3,   -1,  -5,  -9,  -4,  -8,  -4,
            /*C*/ 0, -5,  -1,  -6,  -3,  2,   -2,  -6,  -10, -5,  -9,  -5,
        },
        std::vector{
            //      e,  A,  C,  G,  T,  A,  A,  A,  A,  C,  G,  T
            /*e*/ N, L,   l,   l,   l,   l,   l,   l,   l,   l,   l,   l,
            /*T*/ N, DUL, l,   l,   DUl, l,   l,   l,   l,   l,   l,   DUl,
            /*T*/ N, DUL, DUl, l,   DUl, DUl, l,   l,   l,   l,   l,   DUl,
            /*T*/ N, DUL, DUl, l,   DUl, DUl, l,   l,   l,   l,   l,   DUl,
            /*T*/ N, DUL, DUl, l,   DUl, DUl, l,   l,   l,   l,   l,   DUl,
            /*T*/ N, DUL, DUl, l,   DUl, DUl, l,   l,   l,   l,   l,   DUl,
            /*A*/ N, DUL, L,   l,   l,   DUl, DUl, DUl, DUl, l,   l,   l,
            /*C*/ N, DUL, DUl, L,   l,   l,   l,   l,   l,   DUl, l,   l,
            /*G*/ N, DuL, Ul,  DUl, L,   l,   l,   l,   l,   l,   DUl, l,
            /*T*/ N, DuL, ul,  Ul,  DUL, L,   l,   l,   l,   l,   l,   DUl,
            /*A*/ N, DuL, uL,  ul,  Ul,  DUL, DUL, DUl, DUl, l,   l,   l,
            /*T*/ N, DUL, Dul, uL,  DuL, UL,  DUL, DUL, DUl, DUl, DUl, DUl,
            /*G*/ N, DuL, ul,  Dul, uL,  uL,  DUL, Dul, DuL, Dul, Dul, Ul,
            /*T*/ N, DuL, ul,  ul,  Dul, uL,  DuL, DUl, Dul, DuL, Ul,  Dul,
            /*C*/ N, DuL, Dul, uL,  ul,  uL,  DuL, Dul, DUl, Dul, L,   Ul,
            /*C*/ N, DuL, Dul, uL,  ul,  uL,  DuL, Dul, Dul, DUl, Dul, uL,
            /*C*/ N, DuL, Dul, DuL, ul,  uL,  DuL, Dul, Dul, Dul, DUl, Dul,
            /*C*/ N, DuL, Dul, DuL, ul,  uL,  DuL, Dul, Dul, Dul, Dul, ul,
            /*C*/ N, DuL, DUl, DuL, ul,  ul,  DuL, Dul, Dul, Dul, Dul, ul,
        }};
}();

} // namespace seqan3::test::alignment::fixture::semi_global::affine::unbanded
