// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <vector>

#include "alignment_fixture.hpp"
#include "semi_global_edit_distance_unbanded.hpp"

#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_max_error.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using seqan3::operator""_dna4;

namespace seqan3::test::alignment::fixture::semi_global::edit_distance::max_errors::unbanded
{

using namespace seqan3::test::alignment::fixture::semi_global::edit_distance::unbanded;
using seqan3::detail::column_index_type;
using seqan3::detail::row_index_type;

static auto dna4_01_e255 = []()
{
    return alignment_fixture
    {
        "AACCGGTTAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
        -5,
        "AC---CGGTT",
        "ACGTACG-TA",
        dna4_01.front_coordinate,
        dna4_01.back_coordinate,
        dna4_01.score_vector,
        dna4_01.trace_vector
    };
}();

static auto dna4_01_e5 = []()
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
    /*C*/ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*G*/ 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*T*/ 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
    /*A*/ 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1
    };

    return alignment_fixture
    {
        "AACCGGTTAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{5}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
        -5,
        "AC---CGGTT",
        "ACGTACG-TA",
        dna4_01.front_coordinate,
        dna4_01.back_coordinate,
        dna4_01.score_matrix().mask_matrix(masking_matrix),
        dna4_01.trace_matrix().mask_matrix(masking_matrix)
    };
}();

static auto dna4_01_e2 = []()
{
    std::vector<bool> masking_matrix
    {
    //    e, A, A, C, C, G, G, T, T, A, A, C, C, G, G, T, T
    /*e*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*A*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*C*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*G*/ 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0,
    /*T*/ 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0,
    /*A*/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*C*/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*G*/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*T*/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*A*/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

    return alignment_fixture
    {
        "AACCGGTTAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{2}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
        INF,
        "",
        "",
        seqan3::alignment_coordinate{column_index_type{16u}, row_index_type{9u}},
        seqan3::alignment_coordinate{column_index_type{16u}, row_index_type{9u}},
        dna4_01.score_matrix().mask_matrix(masking_matrix),
        dna4_01.trace_matrix().mask_matrix(masking_matrix)
    };
}();

static auto dna4_01T_e255 = []()
{
    return alignment_fixture
    {
        "ACGTACGTA"_dna4,
        "AACCGGTTAACCGGTT"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
        -8,
        "A-C-G-T-A-C-G-TA",
        "AACCGGTTAACCGGTT",
        dna4_01T.front_coordinate,
        dna4_01T.back_coordinate,
        dna4_01T.score_vector,
        dna4_01T.trace_vector
    };
}();

static auto dna4_02_e255 = []()
{
    return alignment_fixture
    {
        "AACCGGTAAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
        -4,
        "AC---CGGTA",
        "ACGTACG-TA",
        dna4_02.front_coordinate,
        dna4_02.back_coordinate,
        dna4_02.score_vector,
        dna4_02.trace_vector
    };
}();

static auto dna4_02_s10u_15u_e255 = []()
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
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
        -4,
        "AC---CGGTA",
        "ACGTACG-TA",
        dna4_02_s10u_15u.front_coordinate,
        dna4_02_s10u_15u.back_coordinate,
        dna4_02_s10u_15u.score_vector,
        dna4_02_s10u_15u.trace_vector
    };
}();

static auto dna4_02_s10u_15u_e4 = []()
{
    std::vector<bool> masking_matrix
    {
    //    e, A, A, C, C, G, G, T, A, A, A, C, C, G, G
    /*e*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*A*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*C*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*G*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*T*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*A*/ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*C*/ 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*G*/ 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*T*/ 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1,
    /*A*/ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0
    };

    return alignment_fixture
    {
        // score: INF no alignment
        "AACCGGTAAACCGG"_dna4,
        "ACGTACGTA"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{4}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
        -4,
        "AC---CGGTA",
        "ACGTACG-TA",
        dna4_02_s10u_15u.front_coordinate,
        dna4_02_s10u_15u.back_coordinate,
        dna4_02_s10u_15u.score_matrix().mask_matrix(masking_matrix),
        dna4_02_s10u_15u.trace_matrix().mask_matrix(masking_matrix)
    };
}();

static auto dna4_02_s10u_15u_e3 = []()
{
    std::vector<bool> masking_matrix
    {
    //    e, A, A, C, C, G, G, T, A, A, A, C, C, G, G
    /*e*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*A*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*C*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*G*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*T*/ 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*A*/ 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*C*/ 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0,
    /*G*/ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    /*T*/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    /*A*/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

    return alignment_fixture
    {
        // score: INF no alignment
        "AACCGGTAAACCGG"_dna4,
        "ACGTACGTA"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{3}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
        INF,
        "",
        "",
        seqan3::alignment_coordinate{column_index_type{14u}, row_index_type{9u}},
        seqan3::alignment_coordinate{column_index_type{14u}, row_index_type{9u}},
        dna4_02_s10u_15u.score_matrix().mask_matrix(masking_matrix),
        dna4_02_s10u_15u.trace_matrix().mask_matrix(masking_matrix)
    };
}();

static auto dna4_02_s3u_15u_e255 = []()
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
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
        -0,
        "AC",
        "AC",
        dna4_02_s3u_15u.front_coordinate,
        dna4_02_s3u_15u.back_coordinate,
        dna4_02_s3u_15u.score_vector,
        dna4_02_s3u_15u.trace_vector
    };
}();

static auto dna4_02_s3u_15u_e0 = []()
{
    std::vector<bool> masking_matrix
    {
    //    e, A, A, C, C, G, G, T, A, A, A, C, C, G, G
    /*e*/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    /*A*/ 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0,
    /*C*/ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0
    };

    return alignment_fixture
    {
        // score: 0 (0 deletions, 0 insertion)
        // alignment:
        // AACCGGTAAACCGGTT
        //          ||
        // ---------AC-----
        "AACCGGTAAACCGG"_dna4,
        "AC"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{0}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
        -0,
        "AC",
        "AC",
        dna4_02_s3u_15u.front_coordinate,
        dna4_02_s3u_15u.back_coordinate,
        dna4_02_s3u_15u.score_matrix().mask_matrix(masking_matrix),
        dna4_02_s3u_15u.trace_matrix().mask_matrix(masking_matrix)
    };
}();

static auto dna4_02_s1u_15u_e255 = []()
{
    return alignment_fixture
    {
        // score: 0 - empty alignment
        "AACCGGTAAACCGG"_dna4,
        ""_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
        -0,
        "",
        "",
        dna4_02_s1u_15u.front_coordinate,
        dna4_02_s1u_15u.back_coordinate,
        dna4_02_s1u_15u.score_vector,
        dna4_02_s1u_15u.trace_vector
    };
}();

static auto dna4_02_s1u_15u_e0 = []()
{
    return alignment_fixture
    {
        // score: 0 - empty alignment
        "AACCGGTAAACCGG"_dna4,
        ""_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{0}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
        -0,
        "",
        "",
        dna4_02_s1u_15u.front_coordinate,
        dna4_02_s1u_15u.back_coordinate,
        dna4_02_s1u_15u.score_vector,
        dna4_02_s1u_15u.trace_vector
    };
}();

static auto dna4_01T_s17u_1u_e255 = []()
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
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
        -16,
        "----------------",
        "AACCGGTTAACCGGTT",
        dna4_01T_s17u_1u.front_coordinate,
        dna4_01T_s17u_1u.back_coordinate,
        dna4_01T_s17u_1u.score_vector,
        dna4_01T_s17u_1u.trace_vector
    };
}();

static auto dna4_01T_s17u_1u_e5 = []()
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
    /*T*/ 0,
    /*A*/ 0,
    /*A*/ 0,
    /*C*/ 0,
    /*C*/ 0,
    /*G*/ 0,
    /*G*/ 0,
    /*T*/ 0,
    /*T*/ 0,
    };
    return alignment_fixture
    {
        // score: INF - empty alignment
        ""_dna4,
        "AACCGGTTAACCGGTT"_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{5}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
        INF,
        "",
        "",
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{16u}},
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{16u}},
        dna4_01T_s17u_1u.score_matrix().mask_matrix(masking_matrix),
        dna4_01T_s17u_1u.trace_matrix().mask_matrix(masking_matrix)
    };
}();

static auto dna4_03_e255 = []()
{
    return alignment_fixture
    {
        // score: 0
        ""_dna4,
        ""_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
        -0,
        "",
        "",
        dna4_03.front_coordinate,
        dna4_03.back_coordinate,
        dna4_03.score_vector,
        dna4_03.trace_vector
    };
}();

static auto dna4_03_e0 = []()
{
    return alignment_fixture
    {
        // score: 0
        ""_dna4,
        ""_dna4,
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{0}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
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
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
        -5,
        "UW---WRRII",
        "UWRIUWR-IU",
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
        seqan3::align_cfg::edit | seqan3::align_cfg::max_error{255}
                                | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first},
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
