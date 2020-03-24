// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <vector>

#include "alignment_fixture.hpp"

#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using seqan3::operator""_dna4;

namespace seqan3::test::alignment::fixture::global::edit_distance::unbanded
{

using seqan3::detail::column_index_type;
using seqan3::detail::row_index_type;

static auto dna4_01 = []()
{
    return alignment_fixture
    {
        // score: 8 (7 insertions, 1 substitutions)
        // alignment:
        // AACCGGTTAACCGGTT
        // | | | | | | | |
        // A-C-G-T-A-C-G-TA
        "AACCGGTTAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        seqan3::align_cfg::edit,
        -8,
        "AACCGGTTAACCGGTT",
        "A-C-G-T-A-C-G-TA",
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        seqan3::alignment_coordinate{column_index_type{16u}, row_index_type{9u}},
        std::vector
        {
        //     e,  A,  A,  C,  C,  G,  G,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/ -0, -1, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,-14,-15,-16,
        /*A*/ -1, -0, -1, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,-14,-15,
        /*C*/ -2, -1, -1, -1, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,-14,
        /*G*/ -3, -2, -2, -2, -2, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,
        /*T*/ -4, -3, -3, -3, -3, -3, -3, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,
        /*A*/ -5, -4, -3, -4, -4, -4, -4, -4, -4, -4, -5, -6, -7, -8, -9,-10,-11,
        /*C*/ -6, -5, -4, -3, -4, -5, -5, -5, -5, -5, -5, -5, -6, -7, -8, -9,-10,
        /*G*/ -7, -6, -5, -4, -4, -4, -5, -6, -6, -6, -6, -6, -6, -6, -7, -8, -9,
        /*T*/ -8, -7, -6, -5, -5, -5, -5, -5, -6, -7, -7, -7, -7, -7, -7, -7, -8,
        /*A*/ -9, -8, -7, -6, -6, -6, -6, -6, -6, -6, -7, -8, -8, -8, -8, -8, -8
        },
        std::vector
        {
        //      e,  A,  A,  C,  C,  G,  G,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/ N  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*A*/ u  ,D  ,Dl ,l  ,l  ,l  ,l  ,l  ,l  ,Dl ,Dl ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*C*/ u  ,u  ,D  ,D  ,Dl ,l  ,l  ,l  ,l  ,l  ,l  ,Dl ,Dl ,l  ,l  ,l  ,l  ,
        /*G*/ u  ,u  ,Du ,Du ,D  ,D  ,Dl ,l  ,l  ,l  ,l  ,l  ,l  ,Dl ,Dl ,l  ,l  ,
        /*T*/ u  ,u  ,Du ,Du ,Du ,Du ,D  ,D  ,Dl ,l  ,l  ,l  ,l  ,l  ,l  ,Dl ,Dl ,
        /*A*/ u  ,Du ,D  ,Dul,Du ,Du ,Du ,Du ,D  ,D  ,Dl ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*C*/ u  ,u  ,u  ,D  ,Dl ,Dul,Du ,Du ,Du ,Du ,D  ,D  ,Dl ,l  ,l  ,l  ,l  ,
        /*G*/ u  ,u  ,u  ,u  ,D  ,D  ,Dl ,Dul,Du ,Du ,Du ,Du ,D  ,D  ,Dl ,l  ,l  ,
        /*T*/ u  ,u  ,u  ,u  ,Du ,Du ,D  ,D  ,Dl ,Dul,Du ,Du ,Du ,Du ,D  ,D  ,Dl ,
        /*A*/ u  ,Du ,Du ,u  ,Du ,Du ,Du ,Du ,D  ,D  ,Dl ,Dul,Du ,Du ,Du ,Du ,D
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
        seqan3::align_cfg::edit,
        -8,
        "A-C-G-T-A-C-G-TA",
        "AACCGGTTAACCGGTT",
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        seqan3::alignment_coordinate{column_index_type{9u}, row_index_type{16u}},
        dna4_01.score_matrix().transpose_matrix(),
        dna4_01.trace_matrix().transpose_matrix()
    };
}();

static auto dna4_02 = []()
{
    return alignment_fixture
    {
        // score: 8 (7 insertions, 1 substitutions)
        // alignment:
        // AACCGGTAAACCGGTT
        // | | | || | | |
        // A-C-G-TA--C-G-TA
        "AACCGGTAAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        seqan3::align_cfg::edit,
        -8,
        "AACCGGTAAACCGGTT",
        "A-C-G-TA--C-G-TA",
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        seqan3::alignment_coordinate{column_index_type{16u}, row_index_type{9u}},
        std::vector
        {
        //     e,  A,  A,  C,  C,  G,  G,  T,  A,  A,  A,  C,  C,  G,  G,  T,  T,
        /*e*/ -0, -1, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,-14,-15,-16,
        /*A*/ -1, -0, -1, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,-14,-15,
        /*C*/ -2, -1, -1, -1, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,-14,
        /*G*/ -3, -2, -2, -2, -2, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,
        /*T*/ -4, -3, -3, -3, -3, -3, -3, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,
        /*A*/ -5, -4, -3, -4, -4, -4, -4, -4, -3, -4, -5, -6, -7, -8, -9,-10,-11,
        /*C*/ -6, -5, -4, -3, -4, -5, -5, -5, -4, -4, -5, -5, -6, -7, -8, -9,-10,
        /*G*/ -7, -6, -5, -4, -4, -4, -5, -6, -5, -5, -5, -6, -6, -6, -7, -8, -9,
        /*T*/ -8, -7, -6, -5, -5, -5, -5, -5, -6, -6, -6, -6, -7, -7, -7, -7, -8,
        /*A*/ -9, -8, -7, -6, -6, -6, -6, -6, -5, -6, -6, -7, -7, -8, -8, -8, -8
        },
        std::vector
        {
        //      e,  A,  A,  C,  C,  G,  G,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/ N  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*A*/ u  ,D  ,Dl ,l  ,l  ,l  ,l  ,l  ,Dl ,Dl ,Dl ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*C*/ u  ,u  ,D  ,D  ,Dl ,l  ,l  ,l  ,l  ,l  ,l  ,Dl ,Dl ,l  ,l  ,l  ,l  ,
        /*G*/ u  ,u  ,Du ,Du ,D  ,D  ,Dl ,l  ,l  ,l  ,l  ,l  ,l  ,Dl ,Dl ,l  ,l  ,
        /*T*/ u  ,u  ,Du ,Du ,Du ,Du ,D  ,D  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,Dl ,Dl ,
        /*A*/ u  ,Du ,D  ,Dul,Du ,Du ,Du ,Du ,D  ,Dl ,Dl ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*C*/ u  ,u  ,u  ,D  ,Dl ,Dul,Du ,Du ,u  ,D  ,Dl ,D  ,Dl ,l  ,l  ,l  ,l  ,
        /*G*/ u  ,u  ,u  ,u  ,D  ,D  ,Dl ,Dul,u  ,Du ,D  ,Dul,D  ,D  ,Dl ,l  ,l  ,
        /*T*/ u  ,u  ,u  ,u  ,Du ,Du ,D  ,D  ,ul ,Du ,Du ,D  ,Dul,Du ,D  ,D  ,Dl ,
        /*A*/ u  ,Du ,Du ,u  ,Du ,Du ,Du ,Du ,D  ,Dl ,D  ,Dul,D  ,Dul,Du ,Du ,D
        }
    };
}();

static auto dna4_02_s10u_15u = []()
{
    return alignment_fixture
    {
        // score: 8 (7 insertions, 1 substitutions)
        // alignment:
        // AACCGGTAAACCGG-
        // | | | || | ||
        // A-C-G-TA--C-GTA
        "AACCGGTAAACCGG"_dna4,
        "ACGTACGTA"_dna4,
        seqan3::align_cfg::edit,
        -8,
        "AACCGGTAAACCGG-",
        "A-C-G-TA--C-GTA",
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        seqan3::alignment_coordinate{column_index_type{14u}, row_index_type{9u}},
        dna4_02.score_matrix().sub_matrix(10u, 15u),
        dna4_02.trace_matrix().sub_matrix(10u, 15u)
    };
}();

static auto dna4_02_s3u_15u = []()
{
    return alignment_fixture
    {
        // score: 12 (12 insertions)
        // alignment:
        // AACCGGTAAACCGG
        // | |
        // A-C-----------
        "AACCGGTAAACCGG"_dna4,
        "AC"_dna4,
        seqan3::align_cfg::edit,
        -12,
        "AACCGGTAAACCGG",
        "A-C-----------",
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        seqan3::alignment_coordinate{column_index_type{14u}, row_index_type{2u}},
        dna4_02.score_matrix().sub_matrix(3u, 15u),
        dna4_02.trace_matrix().sub_matrix(3u, 15u)
    };
}();

static auto dna4_02_s1u_15u = []()
{
    return alignment_fixture
    {
        // score: 14 (14 deletetions)
        // alignment:
        // AACCGGTAAACCGG
        //
        // --------------
        "AACCGGTAAACCGG"_dna4,
        ""_dna4,
        seqan3::align_cfg::edit,
        -14,
        "AACCGGTAAACCGG",
        "--------------",
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        seqan3::alignment_coordinate{column_index_type{14u}, row_index_type{0u}},
        dna4_02.score_matrix().sub_matrix(1u, 15u),
        dna4_02.trace_matrix().sub_matrix(1u, 15u)
    };
}();

static auto dna4_02T_s15u_1u = []()
{
    return alignment_fixture
    {
        // score: 14 (14 insertions)
        // alignment:
        // --------------
        //
        // AACCGGTAAACCGG
        ""_dna4,
        "AACCGGTAAACCGG"_dna4,
        seqan3::align_cfg::edit,
        -14,
        "--------------",
        "AACCGGTAAACCGG",
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{14u}},
        dna4_02.score_matrix().transpose_matrix().sub_matrix(15u, 1u),
        dna4_02.trace_matrix().transpose_matrix().sub_matrix(15u, 1u)
    };
}();

static auto dna4_03 = []()
{
    return alignment_fixture
    {
        // score: 0
        ""_dna4,
        ""_dna4,
        seqan3::align_cfg::edit,
        -0,
        "",
        "",
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        std::vector<int>{0},
        std::vector<seqan3::detail::trace_directions>{N}
    };
}();

static auto aa27_01 = []()
{
    return alignment_fixture
    {
        // score: 8 (7 insertions, 1 substitutions)
        // alignment:
        // UUWWRRIIUUWWRRII
        // | | | | | | | |
        // U-W-R-I-U-W-R-IU
        "UUWWRRIIUUWWRRII"_aa27,
        "UWRIUWRIU"_aa27,
        seqan3::align_cfg::edit,
        -8,
        "UUWWRRIIUUWWRRII",
        "U-W-R-I-U-W-R-IU",
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        seqan3::alignment_coordinate{column_index_type{16u}, row_index_type{9u}},
        std::vector
        {
        //     e,  U,  U,  W,  W,  R,  R,  I,  I,  U,  U,  W,  W,  R,  R,  I,  I
        /*e*/ -0, -1, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,-14,-15,-16,
        /*U*/ -1, -0, -1, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,-14,-15,
        /*W*/ -2, -1, -1, -1, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,-14,
        /*R*/ -3, -2, -2, -2, -2, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,
        /*I*/ -4, -3, -3, -3, -3, -3, -3, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,
        /*U*/ -5, -4, -3, -4, -4, -4, -4, -4, -4, -4, -5, -6, -7, -8, -9,-10,-11,
        /*W*/ -6, -5, -4, -3, -4, -5, -5, -5, -5, -5, -5, -5, -6, -7, -8, -9,-10,
        /*R*/ -7, -6, -5, -4, -4, -4, -5, -6, -6, -6, -6, -6, -6, -6, -7, -8, -9,
        /*I*/ -8, -7, -6, -5, -5, -5, -5, -5, -6, -7, -7, -7, -7, -7, -7, -7, -8,
        /*U*/ -9, -8, -7, -6, -6, -6, -6, -6, -6, -6, -7, -8, -8, -8, -8, -8, -8
        },
        std::vector
        {
        //      e,  U,  U,  W,  W,  R,  R,  I,  I,  U,  U,  W,  W,  R,  R,  I,  I
        /*e*/ N  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*u*/ u  ,D  ,Dl ,l  ,l  ,l  ,l  ,l  ,l  ,Dl ,Dl ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*W*/ u  ,u  ,D  ,D  ,Dl ,l  ,l  ,l  ,l  ,l  ,l  ,Dl ,Dl ,l  ,l  ,l  ,l  ,
        /*R*/ u  ,u  ,Du ,Du ,D  ,D  ,Dl ,l  ,l  ,l  ,l  ,l  ,l  ,Dl ,Dl ,l  ,l  ,
        /*I*/ u  ,u  ,Du ,Du ,Du ,Du ,D  ,D  ,Dl ,l  ,l  ,l  ,l  ,l  ,l  ,Dl ,Dl ,
        /*u*/ u  ,Du ,D  ,Dul,Du ,Du ,Du ,Du ,D  ,D  ,Dl ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*W*/ u  ,u  ,u  ,D  ,Dl ,Dul,Du ,Du ,Du ,Du ,D  ,D  ,Dl ,l  ,l  ,l  ,l  ,
        /*R*/ u  ,u  ,u  ,u  ,D  ,D  ,Dl ,Dul,Du ,Du ,Du ,Du ,D  ,D  ,Dl ,l  ,l  ,
        /*I*/ u  ,u  ,u  ,u  ,Du ,Du ,D  ,D  ,Dl ,Dul,Du ,Du ,Du ,Du ,D  ,D  ,Dl ,
        /*u*/ u  ,Du ,Du ,u  ,Du ,Du ,Du ,Du ,D  ,D  ,Dl ,Dul,Du ,Du ,Du ,Du ,D
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
        seqan3::align_cfg::edit,
        -8,
        "U-W-R-I-U-W-R-IU",
        "UUWWRRIIUUWWRRII",
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        seqan3::alignment_coordinate{column_index_type{9u}, row_index_type{16u}},
        aa27_01.score_matrix().transpose_matrix(),
        aa27_01.trace_matrix().transpose_matrix()
    };
}();

} // namespace seqan3::test::alignment::fixture::global::edit_distance::unbanded
