// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <vector>

#include "alignment_fixture.hpp"
#include "global_edit_distance_unbanded.hpp"

#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_max_error.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using seqan3::operator""_dna4;

/**
 * NOTE: max_errors is a special case where it will produces the same matrix excepts that it will cutoff all scores from
 * the bottom to the top in the matrix until the score does not exceed the allowed error anymore.
 *
 * Thus we can apply some masking to the matrix and get similiar
 */

namespace seqan3::test::alignment::fixture::global::edit_distance::max_errors::unbanded
{

using namespace seqan3::test::alignment::fixture::global::edit_distance::unbanded;
using seqan3::detail::column_index_type;
using seqan3::detail::row_index_type;

static auto dna4_01_e255 = []()
{
    return alignment_fixture
    {
        "AACCGGTTAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255},
        -8,
        "AACCGGTTAACCGGTT",
        "A-C-G-T-A-C-G-TA",
        dna4_01.front_coordinate,
        dna4_01.back_coordinate,
        dna4_01.score_vector,
        dna4_01.trace_vector
    };
}();

static auto dna4_01T_e255 = []()
{
    return alignment_fixture
    {
        "ACGTACGTA"_dna4,
        "AACCGGTTAACCGGTT"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255},
        -8,
        "A-C-G-T-A-C-G-TA",
        "AACCGGTTAACCGGTT",
        dna4_01T.front_coordinate,
        dna4_01T.back_coordinate,
        dna4_01T.score_vector,
        dna4_01T.trace_vector
    };
}();


static auto dna4_01_e8 = []()
{
    std::vector<bool> masking_matrix
    {
    //    e, A, A, C, C, G, G, T, T, A, A, C, C, G, G, T, T
    /*e*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*A*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*C*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*G*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*T*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*A*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*C*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*G*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*T*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*A*/ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
    };

    return alignment_fixture
    {
        "AACCGGTTAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{8},
        -8,
        "AACCGGTTAACCGGTT",
        "A-C-G-T-A-C-G-TA",
        dna4_01.front_coordinate,
        dna4_01.back_coordinate,
        dna4_01.score_matrix().mask_matrix(masking_matrix),
        dna4_01.trace_matrix().mask_matrix(masking_matrix)
    };
}();

static auto dna4_01_e7 = []()
{
    std::vector<bool> masking_matrix
    {
    //    e, A, A, C, C, G, G, T, T, A, A, C, C, G, G, T, T
    /*e*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
    /*A*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
    /*C*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
    /*G*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
    /*T*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
    /*A*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
    /*C*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
    /*G*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
    /*T*/ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
    /*A*/ 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0
    };

    // score is inf and has no alignment
    return alignment_fixture
    {
        "AACCGGTTAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{7},
        INF,
        "",
        "",
        seqan3::alignment_coordinate{column_index_type{16u}, row_index_type{9u}},
        seqan3::alignment_coordinate{column_index_type{16u}, row_index_type{9u}},
        dna4_01.score_matrix().mask_matrix(masking_matrix),
        dna4_01.trace_matrix().mask_matrix(masking_matrix)
    };
}();

static auto dna4_01_e5 = []()
{
    std::vector<bool> masking_matrix
    {
    //    e, A, A, C, C, G, G, T, T, A, A, C, C, G, G, T, T
    /*e*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
    /*A*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
    /*C*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
    /*G*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
    /*T*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
    /*A*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
    /*C*/ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
    /*G*/ 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*T*/ 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*A*/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

    // score is inf and has no alignment
    return alignment_fixture
    {
        "AACCGGTTAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{5},
        INF,
        "",
        "",
        seqan3::alignment_coordinate{column_index_type{16u}, row_index_type{9u}},
        seqan3::alignment_coordinate{column_index_type{16u}, row_index_type{9u}},
        dna4_01.score_matrix().mask_matrix(masking_matrix),
        dna4_01.trace_matrix().mask_matrix(masking_matrix)
    };
}();

static auto dna4_02_e255 = []()
{
    return alignment_fixture
    {
        "AACCGGTAAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255},
        -8,
        "AACCGGTAAACCGGTT",
        "A-C-G-TA--C-G-TA",
        dna4_02.front_coordinate,
        dna4_02.back_coordinate,
        dna4_02.score_vector,
        dna4_02.trace_vector
    };
}();

static auto dna4_02_e8 = []()
{
    std::vector<bool> masking_matrix
    {
    //    e, A, A, C, C, G, G, T, A, A, A, C, C, G, G, T, T,
    /*e*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*A*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*C*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*G*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*T*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*A*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*C*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*G*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*T*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*A*/ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
    };

    return alignment_fixture
    {
        "AACCGGTAAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{8},
        -8,
        "AACCGGTAAACCGGTT",
        "A-C-G-TA--C-G-TA",
        dna4_02.front_coordinate,
        dna4_02.back_coordinate,
        dna4_02.score_matrix().mask_matrix(masking_matrix),
        dna4_02.trace_matrix().mask_matrix(masking_matrix)
    };
}();

static auto dna4_02_e4 = []()
{
    std::vector<bool> masking_matrix
    {
    //    e, A, A, C, C, G, G, T, A, A, A, C, C, G, G, T, T,
    /*e*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
    /*A*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
    /*C*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
    /*G*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
    /*T*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
    /*A*/ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
    /*C*/ 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
    /*G*/ 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*T*/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*A*/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

    return alignment_fixture
    {
        "AACCGGTAAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{4},
        INF,
        "",
        "",
        seqan3::alignment_coordinate{column_index_type{16u}, row_index_type{9u}},
        seqan3::alignment_coordinate{column_index_type{16u}, row_index_type{9u}},
        dna4_02.score_matrix().mask_matrix(masking_matrix),
        dna4_02.trace_matrix().mask_matrix(masking_matrix)
    };
}();

static auto dna4_02_s10u_15u_e7 = []()
{
    std::vector<bool> masking_matrix
    {
    //    e, A, A, C, C, G, G, T, A, A, A, C, C, G, G
    /*e*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*A*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*C*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*G*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*T*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*A*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*C*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*G*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*T*/ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*A*/ 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0
    };

    // score is inf and has no alignment
    return alignment_fixture
    {
        "AACCGGTAAACCGG"_dna4,
        "ACGTACGTA"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{7},
        INF,
        "",
        "",
        seqan3::alignment_coordinate{column_index_type{14u}, row_index_type{9u}},
        seqan3::alignment_coordinate{column_index_type{14u}, row_index_type{9u}},
        dna4_02_s10u_15u.score_matrix().mask_matrix(masking_matrix),
        dna4_02_s10u_15u.trace_matrix().mask_matrix(masking_matrix)
    };
}();

static auto dna4_02_s1u_15u_e255 = []()
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
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255},
        -14,
        "AACCGGTAAACCGG",
        "--------------",
        dna4_02_s1u_15u.front_coordinate,
        dna4_02_s1u_15u.back_coordinate,
        dna4_02_s1u_15u.score_vector,
        dna4_02_s1u_15u.trace_vector
    };
}();

static auto dna4_02_s1u_15u_e5 = []()
{
    std::vector<bool> masking_matrix
    {
    //    e, A, A, C, C, G, G, T, A, A, A, C, C, G, G
    /*e*/ 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

    return alignment_fixture
    {
        // score is inf and has no alignment
        "AACCGGTAAACCGG"_dna4,
        ""_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{5},
        INF,
        "",
        "",
        seqan3::alignment_coordinate{column_index_type{14u}, row_index_type{0u}},
        seqan3::alignment_coordinate{column_index_type{14u}, row_index_type{0u}},
        dna4_02_s1u_15u.score_matrix().mask_matrix(masking_matrix),
        dna4_02_s1u_15u.trace_matrix().mask_matrix(masking_matrix)
    };
}();

static auto dna4_02T_s15u_1u_e255 = []()
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
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255},
        -14,
        "--------------",
        "AACCGGTAAACCGG",
        dna4_02T_s15u_1u.front_coordinate,
        dna4_02T_s15u_1u.back_coordinate,
        dna4_02T_s15u_1u.score_vector,
        dna4_02T_s15u_1u.trace_vector
    };
}();

static auto dna4_02T_s15u_1u_e5 = []()
{
    std::vector<bool> masking_matrix
    {
    //    e,
    /*e*/ 1,
    /*A*/ 1,
    /*A*/ 1,
    /*C*/ 1,
    /*C*/ 1,
    /*G*/ 1,
    /*G*/ 0,
    /*T*/ 0,
    /*A*/ 0,
    /*A*/ 0,
    /*A*/ 0,
    /*C*/ 0,
    /*C*/ 0,
    /*G*/ 0,
    /*G*/ 0,
    };

    return alignment_fixture
    {
        // score is inf and has no alignment
        ""_dna4,
        "AACCGGTAAACCGG"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{5},
        INF,
        "",
        "",
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{14u}},
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{14u}},
        dna4_02T_s15u_1u.score_matrix().mask_matrix(masking_matrix),
        dna4_02T_s15u_1u.trace_matrix().mask_matrix(masking_matrix)
    };
}();

static auto dna4_03_e255 = []()
{
    return alignment_fixture
    {
        // score: 0
        ""_dna4,
        ""_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255},
        -0,
        "",
        "",
        dna4_03.front_coordinate,
        dna4_03.back_coordinate,
        dna4_03.score_vector,
        dna4_03.trace_vector
    };
}();

static auto aa27_01_e255 = []()
{
    return alignment_fixture
    {
        "UUWWRRIIUUWWRRII"_aa27,
        "UWRIUWRIU"_aa27,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255},
        -8,
        "UUWWRRIIUUWWRRII",
        "U-W-R-I-U-W-R-IU",
        aa27_01.front_coordinate,
        aa27_01.back_coordinate,
        aa27_01.score_vector,
        aa27_01.trace_vector
    };
}();

static auto aa27_01T_e255 = []()
{
    return alignment_fixture
    {
        "UWRIUWRIU"_aa27,
        "UUWWRRIIUUWWRRII"_aa27,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255},
        -8,
        "U-W-R-I-U-W-R-IU",
        "UUWWRRIIUUWWRRII",
        aa27_01T.front_coordinate,
        aa27_01T.back_coordinate,
        aa27_01T.score_vector,
        aa27_01T.trace_vector
    };
}();

} // namespace seqan3::test::alignment::fixture::global::edit_distance::unbanded
