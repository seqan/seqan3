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
        /*e*/NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,
        /*A*/U  ,D  ,D  ,DUL,DU ,DU ,DU ,DU ,DU ,D  ,D  ,DUL,DU ,DU ,DU ,DU ,DU ,
        /*C*/U  ,U  ,DU ,D  ,DL ,DUL,DU ,DU ,DU ,U  ,DU ,D  ,DL ,DUL,DU ,DU ,DU ,
        /*G*/U  ,U  ,DU ,U  ,D  ,D  ,DL ,DUL,DU ,U  ,DU ,U  ,D  ,D  ,DL ,DUL,DU ,
        /*T*/U  ,U  ,DU ,U  ,DU ,DU ,D  ,D  ,DL ,U  ,DU ,U  ,DU ,DU ,D  ,D  ,DL ,
        /*A*/U  ,DU ,D  ,U  ,DU ,DU ,DU ,DU ,D  ,D  ,D  ,U  ,DU ,DU ,DU ,DU ,D  ,
        /*C*/U  ,U  ,U  ,D  ,D  ,DUL,DU ,DU ,DU ,DU ,DU ,D  ,D  ,DUL,DU ,DU ,DU ,
        /*G*/U  ,U  ,U  ,U  ,DU ,D  ,DL ,DUL,DU ,DU ,DU ,U  ,DU ,D  ,DL ,DUL,DU ,
        /*T*/U  ,U  ,U  ,U  ,DU ,U  ,D  ,D  ,DL ,DUL,DU ,U  ,DU ,U  ,D  ,D  ,DL ,
        /*A*/U  ,DU ,DU ,U  ,DU ,U  ,DU ,DU ,D  ,D  ,DL ,U  ,DU ,U  ,DU ,DU ,D
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
        /*e*/NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,
        /*A*/U  ,D  ,DUL,DU ,DU ,D  ,DUL,DU ,DU ,D  ,
        /*A*/U  ,DU ,D  ,DUL,DU ,DU ,D  ,DUL,DU ,DU ,
        /*C*/U  ,U  ,D  ,DL ,DUL,U  ,D  ,DL ,DUL,U  ,
        /*C*/U  ,U  ,DU ,D  ,DL ,U  ,DU ,D  ,DL ,U  ,
        /*G*/U  ,U  ,U  ,D  ,DL ,DUL,U  ,D  ,DL ,DUL,
        /*G*/U  ,U  ,U  ,DU ,D  ,DL ,U  ,DU ,D  ,DL ,
        /*T*/U  ,U  ,U  ,U  ,D  ,DL ,DUL,U  ,D  ,DL ,
        /*T*/U  ,U  ,U  ,U  ,DU ,D  ,DL ,U  ,DU ,D  ,
        /*A*/U  ,DU ,U  ,U  ,U  ,D  ,DL ,DUL,U  ,D  ,
        /*A*/U  ,DU ,U  ,U  ,U  ,DU ,D  ,DL ,U  ,DU ,
        /*C*/U  ,U  ,DU ,U  ,U  ,U  ,D  ,DL ,DUL,U  ,
        /*C*/U  ,U  ,DU ,U  ,U  ,U  ,DU ,D  ,DL ,U  ,
        /*G*/U  ,U  ,U  ,DU ,U  ,U  ,U  ,D  ,DL ,DUL,
        /*G*/U  ,U  ,U  ,DU ,U  ,U  ,U  ,DU ,D  ,DL ,
        /*T*/U  ,U  ,U  ,U  ,DU ,U  ,U  ,U  ,D  ,DL ,
        /*T*/U  ,U  ,U  ,U  ,DU ,U  ,U  ,U  ,DU ,D
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
        /*e*/NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,
        /*A*/U  ,D  ,D  ,DUL,DU ,DU ,DU ,DU ,D  ,D  ,D  ,DUL,DU ,DU ,DU ,DU ,DU ,
        /*C*/U  ,U  ,DU ,D  ,DL ,DUL,DU ,DU ,U  ,DU ,DU ,D  ,DL ,DUL,DU ,DU ,DU ,
        /*G*/U  ,U  ,DU ,U  ,D  ,D  ,DL ,DUL,U  ,DU ,DU ,U  ,D  ,D  ,DL ,DUL,DU ,
        /*T*/U  ,U  ,DU ,U  ,DU ,DU ,D  ,D  ,UL ,DU ,DU ,U  ,DU ,DU ,D  ,D  ,DL ,
        /*A*/U  ,DU ,D  ,U  ,DU ,DU ,DU ,DU ,D  ,DL ,D  ,U  ,DU ,DU ,DU ,DU ,D  ,
        /*C*/U  ,U  ,U  ,D  ,D  ,DUL,DU ,DU ,U  ,D  ,DUL,D  ,D  ,DUL,DU ,DU ,DU ,
        /*G*/U  ,U  ,U  ,U  ,DU ,D  ,DL ,DUL,U  ,DU ,D  ,U  ,DU ,D  ,DL ,DUL,DU ,
        /*T*/U  ,U  ,U  ,U  ,DU ,U  ,D  ,D  ,UL ,DU ,DU ,DU ,DU ,U  ,D  ,D  ,DL ,
        /*A*/U  ,DU ,DU ,U  ,DU ,U  ,DU ,DU ,D  ,DL ,D  ,DUL,DU ,U  ,DU ,DU ,D
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
        0,
        "AC",
        "AC",
        alignment_coordinate{column_index_type{9u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{11u}, row_index_type{2u}},
        dna4_02.score_matrix().sub_matrix(3u, 15u),
        dna4_02.trace_matrix().sub_matrix(3u, 15u)
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
        /*e*/NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,
        /*U*/U  ,D  ,D  ,DUL,DU ,DU ,DU ,DU ,DU ,D  ,D  ,DUL,DU ,DU ,DU ,DU ,DU ,
        /*W*/U  ,U  ,DU ,D  ,DL ,DUL,DU ,DU ,DU ,U  ,DU ,D  ,DL ,DUL,DU ,DU ,DU ,
        /*R*/U  ,U  ,DU ,U  ,D  ,D  ,DL ,DUL,DU ,U  ,DU ,U  ,D  ,D  ,DL ,DUL,DU ,
        /*I*/U  ,U  ,DU ,U  ,DU ,DU ,D  ,D  ,DL ,U  ,DU ,U  ,DU ,DU ,D  ,D  ,DL ,
        /*U*/U  ,DU ,D  ,U  ,DU ,DU ,DU ,DU ,D  ,D  ,D  ,U  ,DU ,DU ,DU ,DU ,D  ,
        /*W*/U  ,U  ,U  ,D  ,D  ,DUL,DU ,DU ,DU ,DU ,DU ,D  ,D  ,DUL,DU ,DU ,DU ,
        /*R*/U  ,U  ,U  ,U  ,DU ,D  ,DL ,DUL,DU ,DU ,DU ,U  ,DU ,D  ,DL ,DUL,DU ,
        /*I*/U  ,U  ,U  ,U  ,DU ,U  ,D  ,D  ,DL ,DUL,DU ,U  ,DU ,U  ,D  ,D  ,DL ,
        /*U*/U  ,DU ,DU ,U  ,DU ,U  ,DU ,DU ,D  ,D  ,DL ,U  ,DU ,U  ,DU ,DU ,D
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
        /*e*/NON,NON,NON,NON,NON,NON,NON,NON,NON,NON,
        /*U*/U  ,D  ,DUL,DU ,DU ,D  ,DUL,DU ,DU ,D  ,
        /*U*/U  ,DU ,D  ,DUL,DU ,DU ,D  ,DUL,DU ,DU ,
        /*W*/U  ,U  ,D  ,DL ,DUL,U  ,D  ,DL ,DUL,U  ,
        /*W*/U  ,U  ,DU ,D  ,DL ,U  ,DU ,D  ,DL ,U  ,
        /*R*/U  ,U  ,U  ,D  ,DL ,DUL,U  ,D  ,DL ,DUL,
        /*R*/U  ,U  ,U  ,DU ,D  ,DL ,U  ,DU ,D  ,DL ,
        /*I*/U  ,U  ,U  ,U  ,D  ,DL ,DUL,U  ,D  ,DL ,
        /*I*/U  ,U  ,U  ,U  ,DU ,D  ,DL ,U  ,DU ,D  ,
        /*U*/U  ,DU ,U  ,U  ,U  ,D  ,DL ,DUL,U  ,D  ,
        /*U*/U  ,DU ,U  ,U  ,U  ,DU ,D  ,DL ,U  ,DU ,
        /*W*/U  ,U  ,DU ,U  ,U  ,U  ,D  ,DL ,DUL,U  ,
        /*W*/U  ,U  ,DU ,U  ,U  ,U  ,DU ,D  ,DL ,U  ,
        /*R*/U  ,U  ,U  ,DU ,U  ,U  ,U  ,D  ,DL ,DUL,
        /*R*/U  ,U  ,U  ,DU ,U  ,U  ,U  ,DU ,D  ,DL ,
        /*I*/U  ,U  ,U  ,U  ,DU ,U  ,U  ,U  ,D  ,DL ,
        /*I*/U  ,U  ,U  ,U  ,DU ,U  ,U  ,U  ,DU ,D
        }
    };
}();

} // namespace seqan3::test::alignment::fixture::global::edit_distance::unbanded
