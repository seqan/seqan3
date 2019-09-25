// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <vector>

#include "alignment_fixture.hpp"

#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

namespace seqan3::test::alignment::fixture::semi_global::edit_distance::unbanded
{

using detail::column_index_type;
using detail::row_index_type;

static auto dna4_01 = []()
{
    return alignment_fixture
    {
        // score: 5 (3 deletions, 1 insertion, 1 substitutions)
        // alignment:
        // AACCGGTTAAC---CGGTT
        //          ||   || ||
        // ---------ACGTACG-TA
        "AACCGGTTAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        align_cfg::edit | align_cfg::aligned_ends{free_ends_first},
        -5,
        "AC---CGGTT",
        "ACGTACG-TA",
        alignment_coordinate{column_index_type{9u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{16u}, row_index_type{9u}},
        std::vector
        {
        //     e,  A,  A,  C,  C,  G,  G,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/ -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0,
        /*A*/ -1, -0, -0, -1, -1, -1, -1, -1, -1, -0, -0, -1, -1, -1, -1, -1, -1,
        /*C*/ -2, -1, -1, -0, -1, -2, -2, -2, -2, -1, -1, -0, -1, -2, -2, -2, -2,
        /*G*/ -3, -2, -2, -1, -1, -1, -2, -3, -3, -2, -2, -1, -1, -1, -2, -3, -3,
        /*T*/ -4, -3, -3, -2, -2, -2, -2, -2, -3, -3, -3, -2, -2, -2, -2, -2, -3,
        /*A*/ -5, -4, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
        /*C*/ -6, -5, -4, -3, -3, -4, -4, -4, -4, -4, -4, -3, -3, -4, -4, -4, -4,
        /*G*/ -7, -6, -5, -4, -4, -3, -4, -5, -5, -5, -5, -4, -4, -3, -4, -5, -5,
        /*T*/ -8, -7, -6, -5, -5, -4, -4, -4, -5, -6, -6, -5, -5, -4, -4, -4, -5,
        /*A*/ -9, -8, -7, -6, -6, -5, -5, -5, -5, -5, -6, -6, -6, -5, -5, -5, -5,
        },
        std::vector
        {
        //     e,  A,  A,  C,  C,  G,  G,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,
        /*A*/u  ,D  ,D  ,Dul,Du ,Du ,Du ,Du ,Du ,D  ,D  ,Dul,Du ,Du ,Du ,Du ,Du ,
        /*C*/u  ,u  ,Du ,D  ,Dl ,Dul,Du ,Du ,Du ,u  ,Du ,D  ,Dl ,Dul,Du ,Du ,Du ,
        /*G*/u  ,u  ,Du ,u  ,D  ,D  ,Dl ,Dul,Du ,u  ,Du ,u  ,D  ,D  ,Dl ,Dul,Du ,
        /*T*/u  ,u  ,Du ,u  ,Du ,Du ,D  ,D  ,Dl ,u  ,Du ,u  ,Du ,Du ,D  ,D  ,Dl ,
        /*A*/u  ,Du ,D  ,u  ,Du ,Du ,Du ,Du ,D  ,D  ,D  ,u  ,Du ,Du ,Du ,Du ,D  ,
        /*C*/u  ,u  ,u  ,D  ,D  ,Dul,Du ,Du ,Du ,Du ,Du ,D  ,D  ,Dul,Du ,Du ,Du ,
        /*G*/u  ,u  ,u  ,u  ,Du ,D  ,Dl ,Dul,Du ,Du ,Du ,u  ,Du ,D  ,Dl ,Dul,Du ,
        /*T*/u  ,u  ,u  ,u  ,Du ,u  ,D  ,D  ,Dl ,Dul,Du ,u  ,Du ,u  ,D  ,D  ,Dl ,
        /*A*/u  ,Du ,Du ,u  ,Du ,u  ,Du ,Du ,D  ,D  ,Dl ,u  ,Du ,u  ,Du ,Du ,D
        }
    };
}();

static auto dna4_01T = []()
{
    return alignment_fixture
    {
        // score: 8 (7 insertions, 1 substitutions)
        // alignment:
        // A-C-G-T-A-C-G-TA
        // | | | | | | | |
        // AACCGGTTAACCGGTT
        "ACGTACGTA"_dna4,
        "AACCGGTTAACCGGTT"_dna4,
        align_cfg::edit | align_cfg::aligned_ends{free_ends_first},
        -8,
        "A-C-G-T-A-C-G-TA",
        "AACCGGTTAACCGGTT",
        alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{9u}, row_index_type{16u}},
        std::vector
        {
        //     e,  A,  C,  G,  T,  A,  C,  G,  T,  A
        /*e*/ -0, -0, -0, -0, -0, -0, -0, -0, -0, -0,
        /*A*/ -1, -0, -1, -1, -1, -0, -1, -1, -1, -0,
        /*A*/ -2, -1, -1, -2, -2, -1, -1, -2, -2, -1,
        /*C*/ -3, -2, -1, -2, -3, -2, -1, -2, -3, -2,
        /*C*/ -4, -3, -2, -2, -3, -3, -2, -2, -3, -3,
        /*G*/ -5, -4, -3, -2, -3, -4, -3, -2, -3, -4,
        /*G*/ -6, -5, -4, -3, -3, -4, -4, -3, -3, -4,
        /*T*/ -7, -6, -5, -4, -3, -4, -5, -4, -3, -4,
        /*T*/ -8, -7, -6, -5, -4, -4, -5, -5, -4, -4,
        /*A*/ -9, -8, -7, -6, -5, -4, -5, -6, -5, -4,
        /*A*/-10, -9, -8, -7, -6, -5, -5, -6, -6, -5,
        /*C*/-11,-10, -9, -8, -7, -6, -5, -6, -7, -6,
        /*C*/-12,-11,-10, -9, -8, -7, -6, -6, -7, -7,
        /*G*/-13,-12,-11,-10, -9, -8, -7, -6, -7, -8,
        /*G*/-14,-13,-12,-11,-10, -9, -8, -7, -7, -8,
        /*T*/-15,-14,-13,-12,-11,-10, -9, -8, -7, -8,
        /*T*/-16,-15,-14,-13,-12,-11,-10, -9, -8, -8
        },
        std::vector
        {
        //     e,  A,  C,  G,  T,  A,  C,  G,  T,  A
        /*e*/N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,
        /*A*/u  ,D  ,Dul,Du ,Du ,D  ,Dul,Du ,Du ,D  ,
        /*A*/u  ,Du ,D  ,Dul,Du ,Du ,D  ,Dul,Du ,Du ,
        /*C*/u  ,u  ,D  ,Dl ,Dul,u  ,D  ,Dl ,Dul,u  ,
        /*C*/u  ,u  ,Du ,D  ,Dl ,u  ,Du ,D  ,Dl ,u  ,
        /*G*/u  ,u  ,u  ,D  ,Dl ,Dul,u  ,D  ,Dl ,Dul,
        /*G*/u  ,u  ,u  ,Du ,D  ,Dl ,u  ,Du ,D  ,Dl ,
        /*T*/u  ,u  ,u  ,u  ,D  ,Dl ,Dul,u  ,D  ,Dl ,
        /*T*/u  ,u  ,u  ,u  ,Du ,D  ,Dl ,u  ,Du ,D  ,
        /*A*/u  ,Du ,u  ,u  ,u  ,D  ,Dl ,Dul,u  ,D  ,
        /*A*/u  ,Du ,u  ,u  ,u  ,Du ,D  ,Dl ,u  ,Du ,
        /*C*/u  ,u  ,Du ,u  ,u  ,u  ,D  ,Dl ,Dul,u  ,
        /*C*/u  ,u  ,Du ,u  ,u  ,u  ,Du ,D  ,Dl ,u  ,
        /*G*/u  ,u  ,u  ,Du ,u  ,u  ,u  ,D  ,Dl ,Dul,
        /*G*/u  ,u  ,u  ,Du ,u  ,u  ,u  ,Du ,D  ,Dl ,
        /*T*/u  ,u  ,u  ,u  ,Du ,u  ,u  ,u  ,D  ,Dl ,
        /*T*/u  ,u  ,u  ,u  ,Du ,u  ,u  ,u  ,Du ,D
        }
    };
}();

static auto dna4_02 = []()
{
    return alignment_fixture
    {
        // score: 4 (3 deletions, 1 insertion)
        // alignment:
        // AAC---CGGTAAAC---CGGTT
        //  ||   || ||
        // -ACGTACG-TA-----------
        "AACCGGTAAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        align_cfg::edit | align_cfg::aligned_ends{free_ends_first},
        -4,
        "AC---CGGTA",
        "ACGTACG-TA",
        alignment_coordinate{column_index_type{1u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{8u}, row_index_type{9u}},
        std::vector
        {
        //     e,  A,  A,  C,  C,  G,  G,  T,  A,  A,  A,  C,  C,  G,  G,  T,  T,
        /*e*/ -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0,
        /*A*/ -1, -0, -0, -1, -1, -1, -1, -1, -0, -0, -0, -1, -1, -1, -1, -1, -1,
        /*C*/ -2, -1, -1, -0, -1, -2, -2, -2, -1, -1, -1, -0, -1, -2, -2, -2, -2,
        /*G*/ -3, -2, -2, -1, -1, -1, -2, -3, -2, -2, -2, -1, -1, -1, -2, -3, -3,
        /*T*/ -4, -3, -3, -2, -2, -2, -2, -2, -3, -3, -3, -2, -2, -2, -2, -2, -3,
        /*A*/ -5, -4, -3, -3, -3, -3, -3, -3, -2, -3, -3, -3, -3, -3, -3, -3, -3,
        /*C*/ -6, -5, -4, -3, -3, -4, -4, -4, -3, -3, -4, -3, -3, -4, -4, -4, -4,
        /*G*/ -7, -6, -5, -4, -4, -3, -4, -5, -4, -4, -4, -4, -4, -3, -4, -5, -5,
        /*T*/ -8, -7, -6, -5, -5, -4, -4, -4, -5, -5, -5, -5, -5, -4, -4, -4, -5,
        /*A*/ -9, -8, -7, -6, -6, -5, -5, -5, -4, -5, -5, -6, -6, -5, -5, -5, -5
        },
        std::vector
        {
        //     e,  A,  A,  C,  C,  G,  G,  T,  A,  A,  A,  C,  C,  G,  G,  T,  T,
        /*e*/N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,
        /*A*/u  ,D  ,D  ,Dul,Du ,Du ,Du ,Du ,D  ,D  ,D  ,Dul,Du ,Du ,Du ,Du ,Du ,
        /*C*/u  ,u  ,Du ,D  ,Dl ,Dul,Du ,Du ,u  ,Du ,Du ,D  ,Dl ,Dul,Du ,Du ,Du ,
        /*G*/u  ,u  ,Du ,u  ,D  ,D  ,Dl ,Dul,u  ,Du ,Du ,u  ,D  ,D  ,Dl ,Dul,Du ,
        /*T*/u  ,u  ,Du ,u  ,Du ,Du ,D  ,D  ,ul ,Du ,Du ,u  ,Du ,Du ,D  ,D  ,Dl ,
        /*A*/u  ,Du ,D  ,u  ,Du ,Du ,Du ,Du ,D  ,Dl ,D  ,u  ,Du ,Du ,Du ,Du ,D  ,
        /*C*/u  ,u  ,u  ,D  ,D  ,Dul,Du ,Du ,u  ,D  ,Dul,D  ,D  ,Dul,Du ,Du ,Du ,
        /*G*/u  ,u  ,u  ,u  ,Du ,D  ,Dl ,Dul,u  ,Du ,D  ,u  ,Du ,D  ,Dl ,Dul,Du ,
        /*T*/u  ,u  ,u  ,u  ,Du ,u  ,D  ,D  ,ul ,Du ,Du ,Du ,Du ,u  ,D  ,D  ,Dl ,
        /*A*/u  ,Du ,Du ,u  ,Du ,u  ,Du ,Du ,D  ,Dl ,D  ,Dul,Du ,u  ,Du ,Du ,D
        }
    };
}();

static auto dna4_02_s10u_15u = []()
{
    return alignment_fixture
    {
        // score: 4 (3 deletions, 1 insertion)
        // alignment:
        // AAC---CGGTAAAC---CGGTT
        //  ||   || ||
        // -ACGTACG-TA-----------
        "AACCGGTAAACCGG"_dna4,
        "ACGTACGTA"_dna4,
        align_cfg::edit | align_cfg::aligned_ends{free_ends_first},
        -4,
        "AC---CGGTA",
        "ACGTACG-TA",
        alignment_coordinate{column_index_type{1u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{8u}, row_index_type{9u}},
        dna4_02.score_matrix().sub_matrix(10u, 15u),
        dna4_02.trace_matrix().sub_matrix(10u, 15u)
    };
}();

static auto dna4_02_s3u_15u = []()
{
    return alignment_fixture
    {
        // score: 0 (0 deletions, 0 insertion)
        // alignment:
        // AACCGGTAAACCGGTT
        //          ||
        // ---------AC-----
        "AACCGGTAAACCGG"_dna4,
        "AC"_dna4,
        align_cfg::edit | align_cfg::aligned_ends{free_ends_first},
        -0,
        "AC",
        "AC",
        alignment_coordinate{column_index_type{9u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{11u}, row_index_type{2u}},
        dna4_02.score_matrix().sub_matrix(3u, 15u),
        dna4_02.trace_matrix().sub_matrix(3u, 15u)
    };
}();

static auto dna4_02_s1u_15u = []()
{
    return alignment_fixture
    {
        // score: 0 - empty alignment
        "AACCGGTAAACCGG"_dna4,
        ""_dna4,
        align_cfg::edit | align_cfg::aligned_ends{free_ends_first},
        -0,
        "",
        "",
        alignment_coordinate{column_index_type{14u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{14u}, row_index_type{0u}},
        dna4_02.score_matrix().sub_matrix(1u, 15u),
        dna4_02.trace_matrix().sub_matrix(1u, 15u)
    };
}();

static auto dna4_01T_s17u_1u = []()
{
    return alignment_fixture
    {
        // score: 16 (16 insertions)
        // alignment:
        // ----------------
        //
        // AACCGGTTAACCGGTT
        ""_dna4,
        "AACCGGTTAACCGGTT"_dna4,
        align_cfg::edit | align_cfg::aligned_ends{free_ends_first},
        -16,
        "----------------",
        "AACCGGTTAACCGGTT",
        alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{0u}, row_index_type{16u}},
        dna4_01T.score_matrix().sub_matrix(17u, 1u),
        dna4_01T.trace_matrix().sub_matrix(17u, 1u)
    };
}();

static auto dna4_03 = []()
{
    return alignment_fixture
    {
        // score: 0
        ""_dna4,
        ""_dna4,
        align_cfg::edit | align_cfg::aligned_ends{free_ends_first},
        -0,
        "",
        "",
        alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        std::vector<int>{0},
        std::vector<detail::trace_directions>{N}
    };
}();

static auto aa27_01 = []()
{
    return alignment_fixture
    {
        // score: 5 (3 deletions, 1 insertion, 1 substitutions)
        // alignment:
        // UUWWRRIIUUW---WRRII
        //          ||   || ||
        // ---------UWRIUWR-IU
        "UUWWRRIIUUWWRRII"_aa27,
        "UWRIUWRIU"_aa27,
        align_cfg::edit | align_cfg::aligned_ends{free_ends_first},
        -5,
        "UW---WRRII",
        "UWRIUWR-IU",
        alignment_coordinate{column_index_type{9u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{16u}, row_index_type{9u}},
        std::vector
        {
        //     e,  U,  U,  W,  W,  R,  R,  I,  I,  U,  U,  W,  W,  R,  R,  I,  I
        /*e*/ -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0, -0,
        /*U*/ -1, -0, -0, -1, -1, -1, -1, -1, -1, -0, -0, -1, -1, -1, -1, -1, -1,
        /*W*/ -2, -1, -1, -0, -1, -2, -2, -2, -2, -1, -1, -0, -1, -2, -2, -2, -2,
        /*R*/ -3, -2, -2, -1, -1, -1, -2, -3, -3, -2, -2, -1, -1, -1, -2, -3, -3,
        /*I*/ -4, -3, -3, -2, -2, -2, -2, -2, -3, -3, -3, -2, -2, -2, -2, -2, -3,
        /*U*/ -5, -4, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3, -3,
        /*W*/ -6, -5, -4, -3, -3, -4, -4, -4, -4, -4, -4, -3, -3, -4, -4, -4, -4,
        /*R*/ -7, -6, -5, -4, -4, -3, -4, -5, -5, -5, -5, -4, -4, -3, -4, -5, -5,
        /*I*/ -8, -7, -6, -5, -5, -4, -4, -4, -5, -6, -6, -5, -5, -4, -4, -4, -5,
        /*U*/ -9, -8, -7, -6, -6, -5, -5, -5, -5, -5, -6, -6, -6, -5, -5, -5, -5,
        },
        std::vector
        {
        //     e,  U,  U,  W,  W,  R,  R,  I,  I,  U,  U,  W,  W,  R,  R,  I,  I
        /*e*/N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,
        /*U*/u  ,D  ,D  ,Dul,Du ,Du ,Du ,Du ,Du ,D  ,D  ,Dul,Du ,Du ,Du ,Du ,Du ,
        /*W*/u  ,u  ,Du ,D  ,Dl ,Dul,Du ,Du ,Du ,u  ,Du ,D  ,Dl ,Dul,Du ,Du ,Du ,
        /*R*/u  ,u  ,Du ,u  ,D  ,D  ,Dl ,Dul,Du ,u  ,Du ,u  ,D  ,D  ,Dl ,Dul,Du ,
        /*I*/u  ,u  ,Du ,u  ,Du ,Du ,D  ,D  ,Dl ,u  ,Du ,u  ,Du ,Du ,D  ,D  ,Dl ,
        /*U*/u  ,Du ,D  ,u  ,Du ,Du ,Du ,Du ,D  ,D  ,D  ,u  ,Du ,Du ,Du ,Du ,D  ,
        /*W*/u  ,u  ,u  ,D  ,D  ,Dul,Du ,Du ,Du ,Du ,Du ,D  ,D  ,Dul,Du ,Du ,Du ,
        /*R*/u  ,u  ,u  ,u  ,Du ,D  ,Dl ,Dul,Du ,Du ,Du ,u  ,Du ,D  ,Dl ,Dul,Du ,
        /*I*/u  ,u  ,u  ,u  ,Du ,u  ,D  ,D  ,Dl ,Dul,Du ,u  ,Du ,u  ,D  ,D  ,Dl ,
        /*U*/u  ,Du ,Du ,u  ,Du ,u  ,Du ,Du ,D  ,D  ,Dl ,u  ,Du ,u  ,Du ,Du ,D
        }
    };
}();

static auto aa27_01T = []()
{
    return alignment_fixture
    {
        // score: 8 (7 insertions, 1 substitutions)
        // alignment:
        // U-W-R-I-U-W-R-IU
        // | | | | | | | |
        // UUWWRRIIUUWWRRII
        "UWRIUWRIU"_aa27,
        "UUWWRRIIUUWWRRII"_aa27,
        align_cfg::edit | align_cfg::aligned_ends{free_ends_first},
        -8,
        "U-W-R-I-U-W-R-IU",
        "UUWWRRIIUUWWRRII",
        alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{9u}, row_index_type{16u}},
        std::vector
        {
        //     e,  U,  W,  R,  I,  U,  W,  R,  I,  U
        /*e*/ -0, -0, -0, -0, -0, -0, -0, -0, -0, -0,
        /*U*/ -1, -0, -1, -1, -1, -0, -1, -1, -1, -0,
        /*U*/ -2, -1, -1, -2, -2, -1, -1, -2, -2, -1,
        /*W*/ -3, -2, -1, -2, -3, -2, -1, -2, -3, -2,
        /*W*/ -4, -3, -2, -2, -3, -3, -2, -2, -3, -3,
        /*R*/ -5, -4, -3, -2, -3, -4, -3, -2, -3, -4,
        /*R*/ -6, -5, -4, -3, -3, -4, -4, -3, -3, -4,
        /*I*/ -7, -6, -5, -4, -3, -4, -5, -4, -3, -4,
        /*I*/ -8, -7, -6, -5, -4, -4, -5, -5, -4, -4,
        /*U*/ -9, -8, -7, -6, -5, -4, -5, -6, -5, -4,
        /*U*/-10, -9, -8, -7, -6, -5, -5, -6, -6, -5,
        /*W*/-11,-10, -9, -8, -7, -6, -5, -6, -7, -6,
        /*W*/-12,-11,-10, -9, -8, -7, -6, -6, -7, -7,
        /*R*/-13,-12,-11,-10, -9, -8, -7, -6, -7, -8,
        /*R*/-14,-13,-12,-11,-10, -9, -8, -7, -7, -8,
        /*I*/-15,-14,-13,-12,-11,-10, -9, -8, -7, -8,
        /*I*/-16,-15,-14,-13,-12,-11,-10, -9, -8, -8
        },
        std::vector
        {
        //     e,  U,  W,  R,  I,  U,  W,  R,  I,  U
        /*e*/N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,
        /*U*/u  ,D  ,Dul,Du ,Du ,D  ,Dul,Du ,Du ,D  ,
        /*U*/u  ,Du ,D  ,Dul,Du ,Du ,D  ,Dul,Du ,Du ,
        /*W*/u  ,u  ,D  ,Dl ,Dul,u  ,D  ,Dl ,Dul,u  ,
        /*W*/u  ,u  ,Du ,D  ,Dl ,u  ,Du ,D  ,Dl ,u  ,
        /*R*/u  ,u  ,u  ,D  ,Dl ,Dul,u  ,D  ,Dl ,Dul,
        /*R*/u  ,u  ,u  ,Du ,D  ,Dl ,u  ,Du ,D  ,Dl ,
        /*I*/u  ,u  ,u  ,u  ,D  ,Dl ,Dul,u  ,D  ,Dl ,
        /*I*/u  ,u  ,u  ,u  ,Du ,D  ,Dl ,u  ,Du ,D  ,
        /*U*/u  ,Du ,u  ,u  ,u  ,D  ,Dl ,Dul,u  ,D  ,
        /*U*/u  ,Du ,u  ,u  ,u  ,Du ,D  ,Dl ,u  ,Du ,
        /*W*/u  ,u  ,Du ,u  ,u  ,u  ,D  ,Dl ,Dul,u  ,
        /*W*/u  ,u  ,Du ,u  ,u  ,u  ,Du ,D  ,Dl ,u  ,
        /*R*/u  ,u  ,u  ,Du ,u  ,u  ,u  ,D  ,Dl ,Dul,
        /*R*/u  ,u  ,u  ,Du ,u  ,u  ,u  ,Du ,D  ,Dl ,
        /*I*/u  ,u  ,u  ,u  ,Du ,u  ,u  ,u  ,D  ,Dl ,
        /*I*/u  ,u  ,u  ,u  ,Du ,u  ,u  ,u  ,Du ,D
        }
    };
}();

} // namespace seqan3::test::alignment::fixture::global::edit_distance::unbanded
