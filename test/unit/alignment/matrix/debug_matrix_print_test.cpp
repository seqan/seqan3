// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/debug_matrix.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3;

namespace seqan3::detail
{
struct debug_matrix_print_test : public ::testing::Test
{
    static constexpr auto inf = std::nullopt;
    std::vector<dna4> sequence1 = "AACACGTTAACCGGTT"_dna4;
    std::vector<dna4> sequence2 = "ACGTACGT"_dna4;

    detail::row_wise_matrix<std::optional<int>> score_matrix
    {
        {
            0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
            1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
            2,  1,  1,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
            3,  2,  2,  2,  2,  3,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13,
            4,  3,  3,  3,  3,  3,  4,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
            5,  4,  3,  4,  3,  4,  4,  4,  4,  4,  5,  6,  7,  8,  9, 10, 11,
            6,  5,  4,  3,  4,  3,  4,  5,  5,  5,  5,  5,  6,  7,  8,  9, 10,
            7,  6,  5,  4,  4,  4,  3,  4,  5,  6,  6,  6,  6,  6,  7,  8,  9,
            inf,7,  6,  5,  5,  5,  4,  3,  4,  5,  6,  7,  7,  7,  7,  7,  8
        }, 9u, 17u
    };

    detail::trace_directions N{},
        D{detail::trace_directions::diagonal},
        L{detail::trace_directions::left},
        U{detail::trace_directions::up},
        DL{D|L}, DU{D|U}, UL{U|L}, DUL{D|U|L};

    detail::row_wise_matrix<detail::trace_directions> trace_matrix
    {
        {
            N,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,
            U,  D,  DL, L,  DL, L,  L,  L,  L,  DL, DL, L,  L,  L,  L,  L,  L,
            U,  U,  D,  D,  L,  DL, L,  L,  L,  L,  L,  DL, DL, L,  L,  L,  L,
            U,  U,  DU, DU, D,  DL, D,  L,  L,  L,  L,  L,  L,  DL, DL, L,  L,
            U,  U,  DU, DU, DU, D,  DUL,D,  DL, L,  L,  L,  L,  L,  L,  DL, DL,
            U,  DU, D,  DUL,D,  DUL,D,  U,  D,  D,  DL, L,  L,  L,  L,  L,  L,
            U,  U,  U,  D,  UL, D,  L,  DUL,DU, DU, D,  D,  DL, L,  L,  L,  L,
            U,  U,  U,  U,  D,  U,  D,  L,  L,  DUL,DU, DU, D,  D,  DL, L,  L,
            N,  U,  U,  U,  DU, DU, U,  D,  DL, L,  L,  DUL,DU, DU, D,  D,  DL
        }, 9u, 17u
    };

    std::string score_matrix_ascii
    {
        " ;  ;  ;  ;  ;  ;  ;  ;  ;  ;  ;  ;  ;  ;  ;  ;  ;  ;\n"
        " ;0 ;1 ;2 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;13;14;15;16;\n"
        " ;1 ;0 ;1 ;2 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;13;14;15;\n"
        " ;2 ;1 ;1 ;1 ;2 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;13;14;\n"
        " ;3 ;2 ;2 ;2 ;2 ;3 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;13;\n"
        " ;4 ;3 ;3 ;3 ;3 ;3 ;4 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;\n"
        " ;5 ;4 ;3 ;4 ;3 ;4 ;4 ;4 ;4 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;\n"
        " ;6 ;5 ;4 ;3 ;4 ;3 ;4 ;5 ;5 ;5 ;5 ;5 ;6 ;7 ;8 ;9 ;10;\n"
        " ;7 ;6 ;5 ;4 ;4 ;4 ;3 ;4 ;5 ;6 ;6 ;6 ;6 ;6 ;7 ;8 ;9 ;\n"
        " ;  ;7 ;6 ;5 ;5 ;5 ;4 ;3 ;4 ;5 ;6 ;7 ;7 ;7 ;7 ;7 ;8 ;\n"
    };

    std::string score_matrix_ascii_with_sequences
    {
        " ;  ;A ;A ;C ;A ;C ;G ;T ;T ;A ;A ;C ;C ;G ;G ;T ;T ;\n"
        " ;0 ;1 ;2 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;13;14;15;16;\n"
        "A;1 ;0 ;1 ;2 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;13;14;15;\n"
        "C;2 ;1 ;1 ;1 ;2 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;13;14;\n"
        "G;3 ;2 ;2 ;2 ;2 ;3 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;13;\n"
        "T;4 ;3 ;3 ;3 ;3 ;3 ;4 ;3 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;12;\n"
        "A;5 ;4 ;3 ;4 ;3 ;4 ;4 ;4 ;4 ;4 ;5 ;6 ;7 ;8 ;9 ;10;11;\n"
        "C;6 ;5 ;4 ;3 ;4 ;3 ;4 ;5 ;5 ;5 ;5 ;5 ;6 ;7 ;8 ;9 ;10;\n"
        "G;7 ;6 ;5 ;4 ;4 ;4 ;3 ;4 ;5 ;6 ;6 ;6 ;6 ;6 ;7 ;8 ;9 ;\n"
        "T;  ;7 ;6 ;5 ;5 ;5 ;4 ;3 ;4 ;5 ;6 ;7 ;7 ;7 ;7 ;7 ;8 ;\n"
    };

    std::string score_matrix_unicode
    {
        " ║ε ║  ║  ║  ║  ║  ║  ║  ║  ║  ║  ║  ║  ║  ║  ║  ║  ║\n"
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
        " ║∞ ║7 ║6 ║5 ║5 ║5 ║4 ║3 ║4 ║5 ║6 ║7 ║7 ║7 ║7 ║7 ║8 ║\n"
    };

    std::string score_matrix_unicode_with_sequences
    {
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
        "T║∞   ║7   ║6   ║5   ║5   ║5   ║4   ║3   ║4   ║5   ║6   ║7   ║7   ║7   ║7   ║7   ║8   ║\n"
    };

    std::string trace_matrix_ascii
    {
        " ;    ;    ;    ;    ;    ;    ;    ;    ;    ;    ;    ;    ;    ;    ;    ;    ;    ;\n"
        " ;N   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;\n"
        " ;U   ;D   ;DL  ;L   ;DL  ;L   ;L   ;L   ;L   ;DL  ;DL  ;L   ;L   ;L   ;L   ;L   ;L   ;\n"
        " ;U   ;U   ;D   ;D   ;L   ;DL  ;L   ;L   ;L   ;L   ;L   ;DL  ;DL  ;L   ;L   ;L   ;L   ;\n"
        " ;U   ;U   ;DU  ;DU  ;D   ;DL  ;D   ;L   ;L   ;L   ;L   ;L   ;L   ;DL  ;DL  ;L   ;L   ;\n"
        " ;U   ;U   ;DU  ;DU  ;DU  ;D   ;DUL ;D   ;DL  ;L   ;L   ;L   ;L   ;L   ;L   ;DL  ;DL  ;\n"
        " ;U   ;DU  ;D   ;DUL ;D   ;DUL ;D   ;U   ;D   ;D   ;DL  ;L   ;L   ;L   ;L   ;L   ;L   ;\n"
        " ;U   ;U   ;U   ;D   ;UL  ;D   ;L   ;DUL ;DU  ;DU  ;D   ;D   ;DL  ;L   ;L   ;L   ;L   ;\n"
        " ;U   ;U   ;U   ;U   ;D   ;U   ;D   ;L   ;L   ;DUL ;DU  ;DU  ;D   ;D   ;DL  ;L   ;L   ;\n"
        " ;N   ;U   ;U   ;U   ;DU  ;DU  ;U   ;D   ;DL  ;L   ;L   ;DUL ;DU  ;DU  ;D   ;D   ;DL  ;\n"
    };

    std::string trace_matrix_ascii_with_sequences
    {
        " ;    ;A   ;A   ;C   ;A   ;C   ;G   ;T   ;T   ;A   ;A   ;C   ;C   ;G   ;G   ;T   ;T   ;\n"
        " ;N   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;\n"
        "A;U   ;D   ;DL  ;L   ;DL  ;L   ;L   ;L   ;L   ;DL  ;DL  ;L   ;L   ;L   ;L   ;L   ;L   ;\n"
        "C;U   ;U   ;D   ;D   ;L   ;DL  ;L   ;L   ;L   ;L   ;L   ;DL  ;DL  ;L   ;L   ;L   ;L   ;\n"
        "G;U   ;U   ;DU  ;DU  ;D   ;DL  ;D   ;L   ;L   ;L   ;L   ;L   ;L   ;DL  ;DL  ;L   ;L   ;\n"
        "T;U   ;U   ;DU  ;DU  ;DU  ;D   ;DUL ;D   ;DL  ;L   ;L   ;L   ;L   ;L   ;L   ;DL  ;DL  ;\n"
        "A;U   ;DU  ;D   ;DUL ;D   ;DUL ;D   ;U   ;D   ;D   ;DL  ;L   ;L   ;L   ;L   ;L   ;L   ;\n"
        "C;U   ;U   ;U   ;D   ;UL  ;D   ;L   ;DUL ;DU  ;DU  ;D   ;D   ;DL  ;L   ;L   ;L   ;L   ;\n"
        "G;U   ;U   ;U   ;U   ;D   ;U   ;D   ;L   ;L   ;DUL ;DU  ;DU  ;D   ;D   ;DL  ;L   ;L   ;\n"
        "T;N   ;U   ;U   ;U   ;DU  ;DU  ;U   ;D   ;DL  ;L   ;L   ;DUL ;DU  ;DU  ;D   ;D   ;DL  ;\n"
    };

    std::string trace_matrix_unicode
    {
        " ║ε  ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║   ║\n"
        " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
        "ε║↺  ║←  ║←  ║←  ║←  ║←  ║←  ║←  ║←  ║←  ║←  ║←  ║←  ║←  ║←  ║←  ║←  ║\n"
        " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
        " ║↑  ║↖  ║↖← ║←  ║↖← ║←  ║←  ║←  ║←  ║↖← ║↖← ║←  ║←  ║←  ║←  ║←  ║←  ║\n"
        " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
        " ║↑  ║↑  ║↖  ║↖  ║←  ║↖← ║←  ║←  ║←  ║←  ║←  ║↖← ║↖← ║←  ║←  ║←  ║←  ║\n"
        " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
        " ║↑  ║↑  ║↖↑ ║↖↑ ║↖  ║↖← ║↖  ║←  ║←  ║←  ║←  ║←  ║←  ║↖← ║↖← ║←  ║←  ║\n"
        " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
        " ║↑  ║↑  ║↖↑ ║↖↑ ║↖↑ ║↖  ║↖↑←║↖  ║↖← ║←  ║←  ║←  ║←  ║←  ║←  ║↖← ║↖← ║\n"
        " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
        " ║↑  ║↖↑ ║↖  ║↖↑←║↖  ║↖↑←║↖  ║↑  ║↖  ║↖  ║↖← ║←  ║←  ║←  ║←  ║←  ║←  ║\n"
        " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
        " ║↑  ║↑  ║↑  ║↖  ║↑← ║↖  ║←  ║↖↑←║↖↑ ║↖↑ ║↖  ║↖  ║↖← ║←  ║←  ║←  ║←  ║\n"
        " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
        " ║↑  ║↑  ║↑  ║↑  ║↖  ║↑  ║↖  ║←  ║←  ║↖↑←║↖↑ ║↖↑ ║↖  ║↖  ║↖← ║←  ║←  ║\n"
        " ╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬═══╬\n"
        " ║↺  ║↑  ║↑  ║↑  ║↖↑ ║↖↑ ║↑  ║↖  ║↖← ║←  ║←  ║↖↑←║↖↑ ║↖↑ ║↖  ║↖  ║↖← ║\n"
    };

    std::string trace_matrix_unicode_with_sequences
    {
        " ║ε   ║A   ║A   ║C   ║A   ║C   ║G   ║T   ║T   ║A   ║A   ║C   ║C   ║G   ║G   ║T   ║T   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "ε║↺   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "A║↑   ║↖   ║↖←  ║←   ║↖←  ║←   ║←   ║←   ║←   ║↖←  ║↖←  ║←   ║←   ║←   ║←   ║←   ║←   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "C║↑   ║↑   ║↖   ║↖   ║←   ║↖←  ║←   ║←   ║←   ║←   ║←   ║↖←  ║↖←  ║←   ║←   ║←   ║←   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "G║↑   ║↑   ║↖↑  ║↖↑  ║↖   ║↖←  ║↖   ║←   ║←   ║←   ║←   ║←   ║←   ║↖←  ║↖←  ║←   ║←   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "T║↑   ║↑   ║↖↑  ║↖↑  ║↖↑  ║↖   ║↖↑← ║↖   ║↖←  ║←   ║←   ║←   ║←   ║←   ║←   ║↖←  ║↖←  ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "A║↑   ║↖↑  ║↖   ║↖↑← ║↖   ║↖↑← ║↖   ║↑   ║↖   ║↖   ║↖←  ║←   ║←   ║←   ║←   ║←   ║←   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "C║↑   ║↑   ║↑   ║↖   ║↑←  ║↖   ║←   ║↖↑← ║↖↑  ║↖↑  ║↖   ║↖   ║↖←  ║←   ║←   ║←   ║←   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "G║↑   ║↑   ║↑   ║↑   ║↖   ║↑   ║↖   ║←   ║←   ║↖↑← ║↖↑  ║↖↑  ║↖   ║↖   ║↖←  ║←   ║←   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "T║↺   ║↑   ║↑   ║↑   ║↖↑  ║↖↑  ║↑   ║↖   ║↖←  ║←   ║←   ║↖↑← ║↖↑  ║↖↑  ║↖   ║↖   ║↖←  ║\n"
    };

    struct debug_matrix_access : public debug_matrix<row_wise_matrix<int>>
    {
        using debug_matrix<row_wise_matrix<int>>::unicode_str_length;
    };

    static constexpr auto unicode_str_length = debug_matrix_access::unicode_str_length;
};

} // namespace seqan3::detail

using typename seqan3::detail::debug_matrix_print_test;
using debug_stream_test = seqan3::detail::debug_matrix_print_test;

TEST_F(debug_matrix_print_test, unicode_str_length)
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
    EXPECT_EQ(unicode_str_length("↑"), 1u);
    EXPECT_EQ(unicode_str_length("↖↑"), 2u);
    EXPECT_EQ(unicode_str_length("←"), 1u);
    EXPECT_EQ(unicode_str_length("↖←"), 2u);
    EXPECT_EQ(unicode_str_length("↑←"), 2u);
    EXPECT_EQ(unicode_str_length("↖↑←"), 3u);
}

TEST_F(debug_matrix_print_test, score_matrix_ascii)
{
    detail::debug_matrix matrix{score_matrix};

    fmtflags2 flags = fmtflags2::default_;

    std::stringstream stream;
    matrix.print(stream, flags);
    EXPECT_EQ(stream.str(), score_matrix_ascii);
}

TEST_F(debug_matrix_print_test, score_matrix_ascii_with_sequences)
{
    detail::debug_matrix matrix{score_matrix, sequence1, sequence2};

    fmtflags2 flags = fmtflags2::default_;
    EXPECT_EQ(matrix.auto_column_width(flags), 2u);

    std::stringstream stream;
    matrix.print(stream, flags);
    EXPECT_EQ(stream.str(), score_matrix_ascii_with_sequences);
}

TEST_F(debug_matrix_print_test, score_matrix_unicode)
{
    detail::debug_matrix matrix{score_matrix};

    fmtflags2 flags = fmtflags2::default_ | fmtflags2::utf8;
    EXPECT_EQ(matrix.auto_column_width(flags), 2u);

    std::stringstream stream;
    matrix.print(stream, flags);
    EXPECT_EQ(stream.str(), score_matrix_unicode);
}

TEST_F(debug_matrix_print_test, score_matrix_unicode_with_sequences)
{
    detail::debug_matrix matrix{score_matrix, sequence1, sequence2};
    matrix.column_width = 4u;

    fmtflags2 flags = fmtflags2::default_ | fmtflags2::utf8;
    EXPECT_EQ(matrix.auto_column_width(flags), 2u);

    std::stringstream stream;
    matrix.print(stream, flags);
    EXPECT_EQ(stream.str(), score_matrix_unicode_with_sequences);
}

TEST_F(debug_matrix_print_test, trace_matrix_ascii)
{
    detail::debug_matrix matrix{trace_matrix};
    matrix.column_width = 4u;

    fmtflags2 flags = fmtflags2::default_;

    EXPECT_EQ(matrix.auto_column_width(flags), 3u);

    std::stringstream stream;
    matrix.print(stream, flags);
    EXPECT_EQ(stream.str(), trace_matrix_ascii);
}

TEST_F(debug_matrix_print_test, trace_matrix_ascii_with_sequences)
{
    detail::debug_matrix matrix{trace_matrix, sequence1, sequence2};
    matrix.column_width = 4u;

    fmtflags2 flags = fmtflags2::default_;

    EXPECT_EQ(matrix.auto_column_width(flags), 3u);

    std::stringstream stream;
    matrix.print(stream, flags);
    EXPECT_EQ(stream.str(), trace_matrix_ascii_with_sequences);
}

TEST_F(debug_matrix_print_test, trace_matrix_unicode)
{
    detail::debug_matrix matrix{trace_matrix};

    fmtflags2 flags = fmtflags2::default_ | fmtflags2::utf8;
    EXPECT_EQ(matrix.auto_column_width(flags), 3u);

    std::stringstream stream;
    matrix.print(stream, flags);
    EXPECT_EQ(stream.str(), trace_matrix_unicode);
}

TEST_F(debug_matrix_print_test, trace_matrix_unicode_with_sequences)
{
    detail::debug_matrix matrix{trace_matrix, sequence1, sequence2};
    matrix.column_width = 4u;

    fmtflags2 flags = fmtflags2::default_ | fmtflags2::utf8;
    EXPECT_EQ(matrix.auto_column_width(flags), 3u);

    std::stringstream stream;
    matrix.print(stream, flags);
    EXPECT_EQ(stream.str(), trace_matrix_unicode_with_sequences);
}

TEST_F(debug_stream_test, score_matrix_ascii)
{
    detail::debug_matrix matrix{score_matrix};

    std::stringstream stream;
    debug_stream_type debug_stream{stream};
    debug_stream << matrix;

    EXPECT_EQ(stream.str(), score_matrix_ascii);
}

TEST_F(debug_stream_test, score_matrix_ascii_with_sequences)
{
    detail::debug_matrix matrix{detail::debug_matrix{score_matrix}, sequence1, sequence2};

    std::stringstream stream;
    debug_stream_type debug_stream{stream};
    debug_stream << matrix;

    EXPECT_EQ(stream.str(), score_matrix_ascii_with_sequences);
}

TEST_F(debug_stream_test, score_matrix_unicode)
{
    detail::debug_matrix matrix{score_matrix};

    std::stringstream stream;
    debug_stream_type debug_stream{stream};
    debug_stream << fmtflags2::utf8 << matrix;

    EXPECT_EQ(stream.str(), score_matrix_unicode);
}

TEST_F(debug_stream_test, score_matrix_unicode_with_sequences)
{
    detail::debug_matrix matrix{detail::debug_matrix{score_matrix}, sequence1, sequence2};
    matrix.column_width = 4u;

    std::stringstream stream;
    debug_stream_type debug_stream{stream};
    debug_stream << fmtflags2::utf8 << matrix;

    EXPECT_EQ(stream.str(), score_matrix_unicode_with_sequences);
}

TEST_F(debug_stream_test, trace_matrix_ascii)
{
    detail::debug_matrix matrix{trace_matrix};
    matrix.column_width = 4u;

    std::stringstream stream;
    debug_stream_type debug_stream{stream};
    debug_stream << matrix;

    EXPECT_EQ(stream.str(), trace_matrix_ascii);
}

TEST_F(debug_stream_test, trace_matrix_ascii_with_sequences)
{
    detail::debug_matrix matrix{detail::debug_matrix{trace_matrix}, sequence1, sequence2};
    matrix.column_width = 4u;

    std::stringstream stream;
    debug_stream_type debug_stream{stream};
    debug_stream << matrix;

    EXPECT_EQ(stream.str(), trace_matrix_ascii_with_sequences);
}

TEST_F(debug_stream_test, trace_matrix_unicode)
{
    detail::debug_matrix matrix{trace_matrix};

    std::stringstream stream;
    debug_stream_type debug_stream{stream};
    debug_stream << fmtflags2::utf8 << matrix;

    EXPECT_EQ(stream.str(), trace_matrix_unicode);
}

TEST_F(debug_stream_test, trace_matrix_unicode_with_sequences)
{
    detail::debug_matrix matrix{detail::debug_matrix{trace_matrix}, sequence1, sequence2};
    matrix.column_width = 4u;

    std::stringstream stream;
    debug_stream_type debug_stream{stream};
    debug_stream << fmtflags2::utf8 << matrix;

    EXPECT_EQ(stream.str(), trace_matrix_unicode_with_sequences);
}
