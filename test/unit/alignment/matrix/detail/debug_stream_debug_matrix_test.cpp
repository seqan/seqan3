// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/detail/debug_matrix.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using seqan3::operator""_dna4;

namespace seqan3::detail
{
struct debug_matrix_stream_test : public ::testing::Test
{
    static constexpr auto inf = std::nullopt;
    std::vector<seqan3::dna4> sequence1 = "AACACGTTAACCGGTT"_dna4;
    std::vector<seqan3::dna4> sequence2 = "ACGTACGT"_dna4;

    seqan3::detail::row_wise_matrix<std::optional<int>> score_matrix{
        seqan3::detail::number_rows{9u},
        seqan3::detail::number_cols{17u},
        std::vector<std::optional<int>>{
            0, 1, 2,  3,  4,  5,  6,   7,  8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 0, 1, 2,  3,  4,  5,  6,  7,
            8, 9, 10, 11, 12, 13, 14,  15, 2, 1, 1,  1,  2,  3,  4,  5,  6,  7, 8, 9, 10, 11, 12, 13, 14, 3,
            2, 2, 2,  2,  3,  3,  4,   5,  6, 7, 8,  9,  10, 11, 12, 13, 4,  3, 3, 3, 3,  3,  4,  3,  4,  5,
            6, 7, 8,  9,  10, 11, 12,  5,  4, 3, 4,  3,  4,  4,  4,  4,  4,  5, 6, 7, 8,  9,  10, 11, 6,  5,
            4, 3, 4,  3,  4,  5,  5,   5,  5, 5, 6,  7,  8,  9,  10, 7,  6,  5, 4, 4, 4,  3,  4,  5,  6,  6,
            6, 6, 6,  7,  8,  9,  inf, 7,  6, 5, 5,  5,  4,  3,  4,  5,  6,  7, 7, 7, 7,  7,  8}};

    seqan3::detail::trace_directions N{}, D{seqan3::detail::trace_directions::diagonal},
        L{seqan3::detail::trace_directions::left}, U{seqan3::detail::trace_directions::up}, DL{D | L}, DU{D | U},
        UL{U | L}, DUL{D | U | L};

    seqan3::detail::row_wise_matrix<detail::trace_directions> trace_matrix{
        seqan3::detail::number_rows{9u},
        seqan3::detail::number_cols{17u},
        std::vector{N,   L,  L,   L, L,  L,  L,  L,  L,   L,  L,  L, L,  L, L, L,   L,  U,  D,   DL, L,  DL,
                    L,   L,  L,   L, DL, DL, L,  L,  L,   L,  L,  L, U,  U, D, D,   L,  DL, L,   L,  L,  L,
                    L,   DL, DL,  L, L,  L,  L,  U,  U,   DU, DU, D, DL, D, L, L,   L,  L,  L,   L,  DL, DL,
                    L,   L,  U,   U, DU, DU, DU, D,  DUL, D,  DL, L, L,  L, L, L,   L,  DL, DL,  U,  DU, D,
                    DUL, D,  DUL, D, U,  D,  D,  DL, L,   L,  L,  L, L,  L, U, U,   U,  D,  UL,  D,  L,  DUL,
                    DU,  DU, D,   D, DL, L,  L,  L,  L,   U,  U,  U, U,  D, U, D,   L,  L,  DUL, DU, DU, D,
                    D,   DL, L,   L, N,  U,  U,  U,  DU,  DU, U,  D, DL, L, L, DUL, DU, DU, D,   D,  DL}};

    std::string score_matrix_ascii{" ;  ;  ;  ;  ;  ;  ;  ;  ;  ;  ;  ;  ;  ;  ;  ;  ;  ;\n"
                                   " ;0 ;1 ;2 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;13;14;15;16;\n"
                                   " ;1 ;0 ;1 ;2 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;13;14;15;\n"
                                   " ;2 ;1 ;1 ;1 ;2 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;13;14;\n"
                                   " ;3 ;2 ;2 ;2 ;2 ;3 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;13;\n"
                                   " ;4 ;3 ;3 ;3 ;3 ;3 ;4 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;\n"
                                   " ;5 ;4 ;3 ;4 ;3 ;4 ;4 ;4 ;4 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;\n"
                                   " ;6 ;5 ;4 ;3 ;4 ;3 ;4 ;5 ;5 ;5 ;5 ;5 ;6 ;7 ;8 ;9 ;10;\n"
                                   " ;7 ;6 ;5 ;4 ;4 ;4 ;3 ;4 ;5 ;6 ;6 ;6 ;6 ;6 ;7 ;8 ;9 ;\n"
                                   " ;  ;7 ;6 ;5 ;5 ;5 ;4 ;3 ;4 ;5 ;6 ;7 ;7 ;7 ;7 ;7 ;8 ;\n"};

    std::string score_matrix_ascii_with_sequences{" ;  ;A ;A ;C ;A ;C ;G ;T ;T ;A ;A ;C ;C ;G ;G ;T ;T ;\n"
                                                  " ;0 ;1 ;2 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;13;14;15;16;\n"
                                                  "A;1 ;0 ;1 ;2 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;13;14;15;\n"
                                                  "C;2 ;1 ;1 ;1 ;2 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;13;14;\n"
                                                  "G;3 ;2 ;2 ;2 ;2 ;3 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;13;\n"
                                                  "T;4 ;3 ;3 ;3 ;3 ;3 ;4 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;\n"
                                                  "A;5 ;4 ;3 ;4 ;3 ;4 ;4 ;4 ;4 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;\n"
                                                  "C;6 ;5 ;4 ;3 ;4 ;3 ;4 ;5 ;5 ;5 ;5 ;5 ;6 ;7 ;8 ;9 ;10;\n"
                                                  "G;7 ;6 ;5 ;4 ;4 ;4 ;3 ;4 ;5 ;6 ;6 ;6 ;6 ;6 ;7 ;8 ;9 ;\n"
                                                  "T;  ;7 ;6 ;5 ;5 ;5 ;4 ;3 ;4 ;5 ;6 ;7 ;7 ;7 ;7 ;7 ;8 ;\n"};

    std::string score_matrix_unicode{" ║ε ║  ║  ║  ║  ║  ║  ║  ║  ║  ║  ║  ║  ║  ║  ║  ║  ║\n"
                                     " ╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬\n"
                                     "ε║0 ║1 ║2 ║3 ║4 ║5 ║6 ║7 ║8 ║9 ║10║11║12║13║14║15║16║\n"
                                     " ╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬\n"
                                     " ║1 ║0 ║1 ║2 ║3 ║4 ║5 ║6 ║7 ║8 ║9 ║10║11║12║13║14║15║\n"
                                     " ╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬\n"
                                     " ║2 ║1 ║1 ║1 ║2 ║3 ║4 ║5 ║6 ║7 ║8 ║9 ║10║11║12║13║14║\n"
                                     " ╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬\n"
                                     " ║3 ║2 ║2 ║2 ║2 ║3 ║3 ║4 ║5 ║6 ║7 ║8 ║9 ║10║11║12║13║\n"
                                     " ╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬\n"
                                     " ║4 ║3 ║3 ║3 ║3 ║3 ║4 ║3 ║4 ║5 ║6 ║7 ║8 ║9 ║10║11║12║\n"
                                     " ╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬\n"
                                     " ║5 ║4 ║3 ║4 ║3 ║4 ║4 ║4 ║4 ║4 ║5 ║6 ║7 ║8 ║9 ║10║11║\n"
                                     " ╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬\n"
                                     " ║6 ║5 ║4 ║3 ║4 ║3 ║4 ║5 ║5 ║5 ║5 ║5 ║6 ║7 ║8 ║9 ║10║\n"
                                     " ╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬\n"
                                     " ║7 ║6 ║5 ║4 ║4 ║4 ║3 ║4 ║5 ║6 ║6 ║6 ║6 ║6 ║7 ║8 ║9 ║\n"
                                     " ╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬══╬\n"
                                     " ║∞ ║7 ║6 ║5 ║5 ║5 ║4 ║3 ║4 ║5 ║6 ║7 ║7 ║7 ║7 ║7 ║8 ║\n"};

    std::string score_matrix_unicode_with_sequences{
        " ║ε   ║A   ║A   ║C   ║A   ║C   ║G   ║T   ║T   ║A   ║A   ║C   ║C   ║G   ║G   ║T   ║T   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "ε║0   ║1   ║2   ║3   ║4   ║5   ║6   ║7   ║8   ║9   ║10  ║11  ║12  ║13  ║14  ║15  ║16  ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "A║1   ║0   ║1   ║2   ║3   ║4   ║5   ║6   ║7   ║8   ║9   ║10  ║11  ║12  ║13  ║14  ║15  ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "C║2   ║1   ║1   ║1   ║2   ║3   ║4   ║5   ║6   ║7   ║8   ║9   ║10  ║11  ║12  ║13  ║14  ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "G║3   ║2   ║2   ║2   ║2   ║3   ║3   ║4   ║5   ║6   ║7   ║8   ║9   ║10  ║11  ║12  ║13  ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "T║4   ║3   ║3   ║3   ║3   ║3   ║4   ║3   ║4   ║5   ║6   ║7   ║8   ║9   ║10  ║11  ║12  ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "A║5   ║4   ║3   ║4   ║3   ║4   ║4   ║4   ║4   ║4   ║5   ║6   ║7   ║8   ║9   ║10  ║11  ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "C║6   ║5   ║4   ║3   ║4   ║3   ║4   ║5   ║5   ║5   ║5   ║5   ║6   ║7   ║8   ║9   ║10  ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "G║7   ║6   ║5   ║4   ║4   ║4   ║3   ║4   ║5   ║6   ║6   ║6   ║6   ║6   ║7   ║8   ║9   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "T║∞   ║7   ║6   ║5   ║5   ║5   ║4   ║3   ║4   ║5   ║6   ║7   ║7   ║7   ║7   ║7   ║8   ║\n"};

    std::string trace_matrix_ascii{
        " ;    ;    ;    ;    ;    ;    ;    ;    ;    ;    ;    ;    ;    ;    ;    ;    ;    ;\n"
        " ;N   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;\n"
        " ;u   ;D   ;Dl  ;l   ;Dl  ;l   ;l   ;l   ;l   ;Dl  ;Dl  ;l   ;l   ;l   ;l   ;l   ;l   ;\n"
        " ;u   ;u   ;D   ;D   ;l   ;Dl  ;l   ;l   ;l   ;l   ;l   ;Dl  ;Dl  ;l   ;l   ;l   ;l   ;\n"
        " ;u   ;u   ;Du  ;Du  ;D   ;Dl  ;D   ;l   ;l   ;l   ;l   ;l   ;l   ;Dl  ;Dl  ;l   ;l   ;\n"
        " ;u   ;u   ;Du  ;Du  ;Du  ;D   ;Dul ;D   ;Dl  ;l   ;l   ;l   ;l   ;l   ;l   ;Dl  ;Dl  ;\n"
        " ;u   ;Du  ;D   ;Dul ;D   ;Dul ;D   ;u   ;D   ;D   ;Dl  ;l   ;l   ;l   ;l   ;l   ;l   ;\n"
        " ;u   ;u   ;u   ;D   ;ul  ;D   ;l   ;Dul ;Du  ;Du  ;D   ;D   ;Dl  ;l   ;l   ;l   ;l   ;\n"
        " ;u   ;u   ;u   ;u   ;D   ;u   ;D   ;l   ;l   ;Dul ;Du  ;Du  ;D   ;D   ;Dl  ;l   ;l   ;\n"
        " ;N   ;u   ;u   ;u   ;Du  ;Du  ;u   ;D   ;Dl  ;l   ;l   ;Dul ;Du  ;Du  ;D   ;D   ;Dl  ;\n"};

    std::string trace_matrix_ascii_with_sequences{
        " ;    ;A   ;A   ;C   ;A   ;C   ;G   ;T   ;T   ;A   ;A   ;C   ;C   ;G   ;G   ;T   ;T   ;\n"
        " ;N   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;l   ;\n"
        "A;u   ;D   ;Dl  ;l   ;Dl  ;l   ;l   ;l   ;l   ;Dl  ;Dl  ;l   ;l   ;l   ;l   ;l   ;l   ;\n"
        "C;u   ;u   ;D   ;D   ;l   ;Dl  ;l   ;l   ;l   ;l   ;l   ;Dl  ;Dl  ;l   ;l   ;l   ;l   ;\n"
        "G;u   ;u   ;Du  ;Du  ;D   ;Dl  ;D   ;l   ;l   ;l   ;l   ;l   ;l   ;Dl  ;Dl  ;l   ;l   ;\n"
        "T;u   ;u   ;Du  ;Du  ;Du  ;D   ;Dul ;D   ;Dl  ;l   ;l   ;l   ;l   ;l   ;l   ;Dl  ;Dl  ;\n"
        "A;u   ;Du  ;D   ;Dul ;D   ;Dul ;D   ;u   ;D   ;D   ;Dl  ;l   ;l   ;l   ;l   ;l   ;l   ;\n"
        "C;u   ;u   ;u   ;D   ;ul  ;D   ;l   ;Dul ;Du  ;Du  ;D   ;D   ;Dl  ;l   ;l   ;l   ;l   ;\n"
        "G;u   ;u   ;u   ;u   ;D   ;u   ;D   ;l   ;l   ;Dul ;Du  ;Du  ;D   ;D   ;Dl  ;l   ;l   ;\n"
        "T;N   ;u   ;u   ;u   ;Du  ;Du  ;u   ;D   ;Dl  ;l   ;l   ;Dul ;Du  ;Du  ;D   ;D   ;Dl  ;\n"};

    std::string trace_matrix_unicode{" ║ε  ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║\n"
                                     " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
                                     "ε║↺  ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║\n"
                                     " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
                                     " ║⇡  ║↖  ║↖⇠ ║⇠  ║↖⇠ ║⇠  ║⇠  ║⇠  ║⇠  ║↖⇠ ║↖⇠ ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║\n"
                                     " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
                                     " ║⇡  ║⇡  ║↖  ║↖  ║⇠  ║↖⇠ ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║↖⇠ ║↖⇠ ║⇠  ║⇠  ║⇠  ║⇠  ║\n"
                                     " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
                                     " ║⇡  ║⇡  ║↖⇡ ║↖⇡ ║↖  ║↖⇠ ║↖  ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║↖⇠ ║↖⇠ ║⇠  ║⇠  ║\n"
                                     " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
                                     " ║⇡  ║⇡  ║↖⇡ ║↖⇡ ║↖⇡ ║↖  ║↖⇡⇠║↖  ║↖⇠ ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║↖⇠ ║↖⇠ ║\n"
                                     " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
                                     " ║⇡  ║↖⇡ ║↖  ║↖⇡⇠║↖  ║↖⇡⇠║↖  ║⇡  ║↖  ║↖  ║↖⇠ ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║⇠  ║\n"
                                     " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
                                     " ║⇡  ║⇡  ║⇡  ║↖  ║⇡⇠ ║↖  ║⇠  ║↖⇡⇠║↖⇡ ║↖⇡ ║↖  ║↖  ║↖⇠ ║⇠  ║⇠  ║⇠  ║⇠  ║\n"
                                     " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
                                     " ║⇡  ║⇡  ║⇡  ║⇡  ║↖  ║⇡  ║↖  ║⇠  ║⇠  ║↖⇡⇠║↖⇡ ║↖⇡ ║↖  ║↖  ║↖⇠ ║⇠  ║⇠  ║\n"
                                     " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
                                     " ║↺  ║⇡  ║⇡  ║⇡  ║↖⇡ ║↖⇡ ║⇡  ║↖  ║↖⇠ ║⇠  ║⇠  ║↖⇡⇠║↖⇡ ║↖⇡ ║↖  ║↖  ║↖⇠ ║\n"};

    std::string trace_matrix_unicode_with_sequences{
        " ║ε   ║A   ║A   ║C   ║A   ║C   ║G   ║T   ║T   ║A   ║A   ║C   ║C   ║G   ║G   ║T   ║T   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "ε║↺   ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "A║⇡   ║↖   ║↖⇠  ║⇠   ║↖⇠  ║⇠   ║⇠   ║⇠   ║⇠   ║↖⇠  ║↖⇠  ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "C║⇡   ║⇡   ║↖   ║↖   ║⇠   ║↖⇠  ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║↖⇠  ║↖⇠  ║⇠   ║⇠   ║⇠   ║⇠   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "G║⇡   ║⇡   ║↖⇡  ║↖⇡  ║↖   ║↖⇠  ║↖   ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║↖⇠  ║↖⇠  ║⇠   ║⇠   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "T║⇡   ║⇡   ║↖⇡  ║↖⇡  ║↖⇡  ║↖   ║↖⇡⇠ ║↖   ║↖⇠  ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║↖⇠  ║↖⇠  ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "A║⇡   ║↖⇡  ║↖   ║↖⇡⇠ ║↖   ║↖⇡⇠ ║↖   ║⇡   ║↖   ║↖   ║↖⇠  ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║⇠   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "C║⇡   ║⇡   ║⇡   ║↖   ║⇡⇠  ║↖   ║⇠   ║↖⇡⇠ ║↖⇡  ║↖⇡  ║↖   ║↖   ║↖⇠  ║⇠   ║⇠   ║⇠   ║⇠   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "G║⇡   ║⇡   ║⇡   ║⇡   ║↖   ║⇡   ║↖   ║⇠   ║⇠   ║↖⇡⇠ ║↖⇡  ║↖⇡  ║↖   ║↖   ║↖⇠  ║⇠   ║⇠   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "T║↺   ║⇡   ║⇡   ║⇡   ║↖⇡  ║↖⇡  ║⇡   ║↖   ║↖⇠  ║⇠   ║⇠   ║↖⇡⇠ ║↖⇡  ║↖⇡  ║↖   ║↖   ║↖⇠  ║\n"};

    struct debug_matrix_access : public debug_matrix<row_wise_matrix<int>>
    {
        using debug_matrix<row_wise_matrix<int>>::unicode_str_length;
    };

    static constexpr auto unicode_str_length = debug_matrix_access::unicode_str_length;
};

} // namespace seqan3::detail

using typename seqan3::detail::debug_matrix_stream_test;
using debug_stream_test = seqan3::detail::debug_matrix_stream_test;
using debug_matrix_printer_test = seqan3::detail::debug_matrix_stream_test;

TEST_F(debug_matrix_stream_test, unicode_str_length)
{
    EXPECT_EQ(unicode_str_length(" "), 1u);
    EXPECT_EQ(unicode_str_length(";"), 1u);
    EXPECT_EQ(unicode_str_length(""), 0u);
    EXPECT_EQ(unicode_str_length("N"), 1u);
    EXPECT_EQ(unicode_str_length("D"), 1u);
    EXPECT_EQ(unicode_str_length("U"), 1u);
    EXPECT_EQ(unicode_str_length("DU"), 2u);
    EXPECT_EQ(unicode_str_length("L"), 1u);
    EXPECT_EQ(unicode_str_length("DL"), 2u);
    EXPECT_EQ(unicode_str_length("UL"), 2u);
    EXPECT_EQ(unicode_str_length("DUL"), 3u);
    EXPECT_EQ(unicode_str_length("|"), 1u);
    EXPECT_EQ(unicode_str_length("-"), 1u);
    EXPECT_EQ(unicode_str_length("/"), 1u);
    EXPECT_EQ(unicode_str_length("INF"), 3u);

    EXPECT_EQ(unicode_str_length("ε"), 1u);
    EXPECT_EQ(unicode_str_length("║"), 1u);
    EXPECT_EQ(unicode_str_length("═"), 1u);
    EXPECT_EQ(unicode_str_length("╬"), 1u);
    EXPECT_EQ(unicode_str_length("∞"), 1u);

    EXPECT_EQ(unicode_str_length("█"), 1u);
    EXPECT_EQ(unicode_str_length("▘"), 1u);
    EXPECT_EQ(unicode_str_length("▝"), 1u);
    EXPECT_EQ(unicode_str_length("▀"), 1u);
    EXPECT_EQ(unicode_str_length("▖"), 1u);
    EXPECT_EQ(unicode_str_length("▌"), 1u);
    EXPECT_EQ(unicode_str_length("▞"), 1u);
    EXPECT_EQ(unicode_str_length("▛"), 1u);
    EXPECT_EQ(unicode_str_length("∞"), 1u);

    EXPECT_EQ(unicode_str_length("⠀"), 1u);
    EXPECT_EQ(unicode_str_length("⠁"), 1u);
    EXPECT_EQ(unicode_str_length("⠈"), 1u);
    EXPECT_EQ(unicode_str_length("⠉"), 1u);
    EXPECT_EQ(unicode_str_length("⠄"), 1u);
    EXPECT_EQ(unicode_str_length("⠅"), 1u);
    EXPECT_EQ(unicode_str_length("⠌"), 1u);
    EXPECT_EQ(unicode_str_length("⠍"), 1u);

    EXPECT_EQ(unicode_str_length("↺"), 1u);
    EXPECT_EQ(unicode_str_length("↖"), 1u);
    EXPECT_EQ(unicode_str_length("⇡"), 1u);
    EXPECT_EQ(unicode_str_length("↖⇡"), 2u);
    EXPECT_EQ(unicode_str_length("⇠"), 1u);
    EXPECT_EQ(unicode_str_length("↖⇠"), 2u);
    EXPECT_EQ(unicode_str_length("⇡⇠"), 2u);
    EXPECT_EQ(unicode_str_length("↖⇡⇠"), 3u);
}

TEST_F(debug_matrix_stream_test, score_matrix_ascii)
{
    seqan3::detail::debug_matrix matrix{score_matrix};

    seqan3::fmtflags2 flags = seqan3::fmtflags2::default_;

    std::stringstream stream;
    matrix.stream_matrix(stream, flags);
    EXPECT_EQ(stream.str(), score_matrix_ascii);
}

TEST_F(debug_matrix_stream_test, score_matrix_ascii_with_sequences)
{
    seqan3::detail::debug_matrix matrix{score_matrix, sequence1, sequence2};

    seqan3::fmtflags2 flags = seqan3::fmtflags2::default_;
    EXPECT_EQ(matrix.auto_column_width(flags), 2u);

    std::stringstream stream;
    matrix.stream_matrix(stream, flags);
    EXPECT_EQ(stream.str(), score_matrix_ascii_with_sequences);
}

TEST_F(debug_matrix_stream_test, score_matrix_unicode)
{
    seqan3::detail::debug_matrix matrix{score_matrix};

    seqan3::fmtflags2 flags = seqan3::fmtflags2::default_ | seqan3::fmtflags2::utf8;
    EXPECT_EQ(matrix.auto_column_width(flags), 2u);

    std::stringstream stream;
    matrix.stream_matrix(stream, flags);
    EXPECT_EQ(stream.str(), score_matrix_unicode);
}

TEST_F(debug_matrix_stream_test, score_matrix_unicode_with_sequences)
{
    seqan3::detail::debug_matrix matrix{score_matrix, sequence1, sequence2};
    matrix.column_width = 4u;

    seqan3::fmtflags2 flags = seqan3::fmtflags2::default_ | seqan3::fmtflags2::utf8;
    EXPECT_EQ(matrix.auto_column_width(flags), 2u);

    std::stringstream stream;
    matrix.stream_matrix(stream, flags);
    EXPECT_EQ(stream.str(), score_matrix_unicode_with_sequences);
}

TEST_F(debug_matrix_stream_test, trace_matrix_ascii)
{
    seqan3::detail::debug_matrix matrix{trace_matrix};
    matrix.column_width = 4u;

    seqan3::fmtflags2 flags = seqan3::fmtflags2::default_;

    EXPECT_EQ(matrix.auto_column_width(flags), 3u);

    std::stringstream stream;
    matrix.stream_matrix(stream, flags);
    EXPECT_EQ(stream.str(), trace_matrix_ascii);
}

TEST_F(debug_matrix_stream_test, trace_matrix_ascii_with_sequences)
{
    seqan3::detail::debug_matrix matrix{trace_matrix, sequence1, sequence2};
    matrix.column_width = 4u;

    seqan3::fmtflags2 flags = seqan3::fmtflags2::default_;

    EXPECT_EQ(matrix.auto_column_width(flags), 3u);

    std::stringstream stream;
    matrix.stream_matrix(stream, flags);
    EXPECT_EQ(stream.str(), trace_matrix_ascii_with_sequences);
}

TEST_F(debug_matrix_stream_test, trace_matrix_unicode)
{
    seqan3::detail::debug_matrix matrix{trace_matrix};

    seqan3::fmtflags2 flags = seqan3::fmtflags2::default_ | seqan3::fmtflags2::utf8;
    EXPECT_EQ(matrix.auto_column_width(flags), 3u);

    std::stringstream stream;
    matrix.stream_matrix(stream, flags);
    EXPECT_EQ(stream.str(), trace_matrix_unicode);
}

TEST_F(debug_matrix_stream_test, trace_matrix_unicode_with_sequences)
{
    seqan3::detail::debug_matrix matrix{trace_matrix, sequence1, sequence2};
    matrix.column_width = 4u;

    seqan3::fmtflags2 flags = seqan3::fmtflags2::default_ | seqan3::fmtflags2::utf8;
    EXPECT_EQ(matrix.auto_column_width(flags), 3u);

    std::stringstream stream;
    matrix.stream_matrix(stream, flags);
    EXPECT_EQ(stream.str(), trace_matrix_unicode_with_sequences);
}

TEST_F(debug_stream_test, score_matrix_ascii)
{
    seqan3::detail::debug_matrix matrix{score_matrix};

    std::stringstream stream;
    seqan3::debug_stream_type debug_stream{stream};
    debug_stream << matrix;

    EXPECT_EQ(stream.str(), score_matrix_ascii);
}

TEST_F(debug_stream_test, score_matrix_ascii_with_sequences)
{
    seqan3::detail::debug_matrix matrix{seqan3::detail::debug_matrix{score_matrix}, sequence1, sequence2};

    std::stringstream stream;
    seqan3::debug_stream_type debug_stream{stream};
    debug_stream << matrix;

    EXPECT_EQ(stream.str(), score_matrix_ascii_with_sequences);
}

TEST_F(debug_stream_test, score_matrix_unicode)
{
    seqan3::detail::debug_matrix matrix{score_matrix};

    std::stringstream stream;
    seqan3::debug_stream_type debug_stream{stream};
    debug_stream << seqan3::fmtflags2::utf8 << matrix;

    EXPECT_EQ(stream.str(), score_matrix_unicode);
}

TEST_F(debug_stream_test, score_matrix_unicode_with_sequences)
{
    seqan3::detail::debug_matrix matrix{seqan3::detail::debug_matrix{score_matrix}, sequence1, sequence2};
    matrix.column_width = 4u;

    std::stringstream stream;
    seqan3::debug_stream_type debug_stream{stream};
    debug_stream << seqan3::fmtflags2::utf8 << matrix;

    EXPECT_EQ(stream.str(), score_matrix_unicode_with_sequences);
}

TEST_F(debug_stream_test, trace_matrix_ascii)
{
    seqan3::detail::debug_matrix matrix{trace_matrix};
    matrix.column_width = 4u;

    std::stringstream stream;
    seqan3::debug_stream_type debug_stream{stream};
    debug_stream << matrix;

    EXPECT_EQ(stream.str(), trace_matrix_ascii);
}

TEST_F(debug_stream_test, trace_matrix_ascii_with_sequences)
{
    seqan3::detail::debug_matrix matrix{seqan3::detail::debug_matrix{trace_matrix}, sequence1, sequence2};
    matrix.column_width = 4u;

    std::stringstream stream;
    seqan3::debug_stream_type debug_stream{stream};
    debug_stream << matrix;

    EXPECT_EQ(stream.str(), trace_matrix_ascii_with_sequences);
}

TEST_F(debug_stream_test, trace_matrix_unicode)
{
    seqan3::detail::debug_matrix matrix{trace_matrix};

    std::stringstream stream;
    seqan3::debug_stream_type debug_stream{stream};
    debug_stream << seqan3::fmtflags2::utf8 << matrix;

    EXPECT_EQ(stream.str(), trace_matrix_unicode);
}

TEST_F(debug_stream_test, trace_matrix_unicode_with_sequences)
{
    seqan3::detail::debug_matrix matrix{seqan3::detail::debug_matrix{trace_matrix}, sequence1, sequence2};
    matrix.column_width = 4u;

    std::stringstream stream;
    seqan3::debug_stream_type debug_stream{stream};
    debug_stream << seqan3::fmtflags2::utf8 << matrix;

    EXPECT_EQ(stream.str(), trace_matrix_unicode_with_sequences);
}

TEST_F(debug_matrix_printer_test, score_matrix_ascii_std_stream)
{
    seqan3::detail::debug_matrix matrix{score_matrix};
    seqan3::alignment_matrix_printer<decltype(matrix)> printer;

    std::stringstream stream;
    printer(stream, matrix);

    EXPECT_EQ(stream.str(), score_matrix_ascii);
}

TEST_F(debug_matrix_printer_test, score_matrix_unicode_debug_stream)
{
    seqan3::detail::debug_matrix matrix{score_matrix};
    seqan3::alignment_matrix_printer<decltype(matrix)> printer;

    std::stringstream stream;
    seqan3::debug_stream_type debug_stream{stream};
    printer(debug_stream << seqan3::fmtflags2::utf8, matrix);

    EXPECT_EQ(stream.str(), score_matrix_unicode);
}
