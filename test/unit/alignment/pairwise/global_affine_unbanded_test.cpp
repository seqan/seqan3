// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>


#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>

#include <seqan3/range/view/to_char.hpp>

#include "fixture/global_affine_unbanded.hpp"

using namespace seqan3;
using namespace seqan3::detail;
using namespace seqan3::fixture;

template <auto _fixture>
struct param : public ::testing::Test
{
    auto fixture() -> decltype(alignment_fixture{*_fixture}) const &
    {
        return *_fixture;
    }
};

template <typename param_t>
class global_affine_unbanded : public param_t
{};

TYPED_TEST_CASE_P(global_affine_unbanded);

using global_affine_unbanded_types
    = ::testing::Types<
        param<&global::affine::unbanded::dna4_01>,
        param<&global::affine::unbanded::dna4_02>//,
        // param<&global::affine::unbanded::dna4_01T>,
        // param<&global::affine::unbanded::dna4_02>,
        // param<&global::affine::unbanded::aa27_01>,
        // param<&global::affine::unbanded::aa27_01T>
    >;

// using semi_global_edit_distance_unbanded_types
//     = ::testing::Types<
//         param<&semi_global::edit_distance::unbanded::dna4_01, uint8_t>,
//         param<&semi_global::edit_distance::unbanded::dna4_01, uint16_t>,
//         param<&semi_global::edit_distance::unbanded::dna4_01, uint32_t>,
//         param<&semi_global::edit_distance::unbanded::dna4_01, uint64_t>,
//
//         param<&semi_global::edit_distance::unbanded::dna4_01T, uint8_t>,
//         param<&semi_global::edit_distance::unbanded::dna4_01T, uint16_t>,
//         param<&semi_global::edit_distance::unbanded::dna4_01T, uint32_t>,
//         param<&semi_global::edit_distance::unbanded::dna4_01T, uint64_t>,
//
//         param<&semi_global::edit_distance::unbanded::dna4_02, uint8_t>,
//         param<&semi_global::edit_distance::unbanded::dna4_02, uint16_t>,
//         param<&semi_global::edit_distance::unbanded::dna4_02, uint32_t>,
//         param<&semi_global::edit_distance::unbanded::dna4_02, uint64_t>,
//
//         param<&semi_global::edit_distance::unbanded::aa27_01, uint8_t>,
//         param<&semi_global::edit_distance::unbanded::aa27_01, uint16_t>,
//         param<&semi_global::edit_distance::unbanded::aa27_01, uint32_t>,
//         param<&semi_global::edit_distance::unbanded::aa27_01, uint64_t>,
//
//         param<&semi_global::edit_distance::unbanded::aa27_01T, uint8_t>,
//         param<&semi_global::edit_distance::unbanded::aa27_01T, uint16_t>,
//         param<&semi_global::edit_distance::unbanded::aa27_01T, uint32_t>,
//         param<&semi_global::edit_distance::unbanded::aa27_01T, uint64_t>
//     >;
//
// using global_edit_distance_max_errors_unbanded_types
//     = ::testing::Types<
//         param<&global::edit_distance::max_errors::unbanded::dna4_01_e255, uint8_t>,
//         param<&global::edit_distance::max_errors::unbanded::dna4_01_e255, uint16_t>,
//         param<&global::edit_distance::max_errors::unbanded::dna4_01_e255, uint32_t>,
//         param<&global::edit_distance::max_errors::unbanded::dna4_01_e255, uint64_t>,
//
//         param<&global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint8_t>,
//         param<&global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint16_t>,
//         param<&global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint32_t>,
//         param<&global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint64_t>,
//
//         param<&global::edit_distance::max_errors::unbanded::dna4_02_e255, uint8_t>,
//         param<&global::edit_distance::max_errors::unbanded::dna4_02_e255, uint16_t>,
//         param<&global::edit_distance::max_errors::unbanded::dna4_02_e255, uint32_t>,
//         param<&global::edit_distance::max_errors::unbanded::dna4_02_e255, uint64_t>,
//
//         param<&global::edit_distance::max_errors::unbanded::aa27_01_e255, uint8_t>,
//         param<&global::edit_distance::max_errors::unbanded::aa27_01_e255, uint16_t>,
//         param<&global::edit_distance::max_errors::unbanded::aa27_01_e255, uint32_t>,
//         param<&global::edit_distance::max_errors::unbanded::aa27_01_e255, uint64_t>,
//
//         param<&global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint8_t>,
//         param<&global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint16_t>,
//         param<&global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint32_t>,
//         param<&global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint64_t>
//     >;
//
// using semi_global_edit_distance_max_errors_unbanded_types
//     = ::testing::Types<
//         param<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e255, uint8_t>,
//         param<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e255, uint16_t>,
//         param<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e255, uint32_t>,
//         param<&semi_global::edit_distance::max_errors::unbanded::dna4_01_e255, uint64_t>,
//
//         param<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint8_t>,
//         param<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint16_t>,
//         param<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint32_t>,
//         param<&semi_global::edit_distance::max_errors::unbanded::dna4_01T_e255, uint64_t>,
//
//         param<&semi_global::edit_distance::max_errors::unbanded::dna4_02_e255, uint8_t>,
//         param<&semi_global::edit_distance::max_errors::unbanded::dna4_02_e255, uint16_t>,
//         param<&semi_global::edit_distance::max_errors::unbanded::dna4_02_e255, uint32_t>,
//         param<&semi_global::edit_distance::max_errors::unbanded::dna4_02_e255, uint64_t>,
//
//         param<&semi_global::edit_distance::max_errors::unbanded::aa27_01_e255, uint8_t>,
//         param<&semi_global::edit_distance::max_errors::unbanded::aa27_01_e255, uint16_t>,
//         param<&semi_global::edit_distance::max_errors::unbanded::aa27_01_e255, uint32_t>,
//         param<&semi_global::edit_distance::max_errors::unbanded::aa27_01_e255, uint64_t>,
//
//         param<&semi_global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint8_t>,
//         param<&semi_global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint16_t>,
//         param<&semi_global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint32_t>,
//         param<&semi_global::edit_distance::max_errors::unbanded::aa27_01T_e255, uint64_t>
//     >;

TYPED_TEST_P(global_affine_unbanded, score)
{
    auto const & fixture = this->fixture();
    // We only compute the score.
    auto align_cfg = fixture.config | align_cfg::result{align_cfg::with_score};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    // TODO Make test work with ranges.
    auto alignment = align_pairwise(std::pair{database, query}, align_cfg);

    EXPECT_EQ((*std::ranges::begin(alignment)).get_score(), fixture.score);
}

// TYPED_TEST_P(global_affine_unbanded, score_matrix)
// {
//     using word_type = typename TypeParam::word_type;
//     auto const & fixture = this->fixture();
//     auto align_cfg = fixture.config;
//
//     std::vector database = fixture.sequence1;
//     std::vector query = fixture.sequence2;
//
//     auto alignment = edit_distance<word_type>(database, query, align_cfg);
//     auto score_matrix = alignment.score_matrix();
//
//     EXPECT_EQ(score_matrix.cols(), database.size()+1);
//     EXPECT_EQ(score_matrix.rows(), query.size()+1);
//     EXPECT_EQ(score_matrix, fixture.score_matrix);
//     EXPECT_EQ(alignment.score(), fixture.score);
// }
//
// TYPED_TEST_P(global_affine_unbanded, trace_matrix)
// {
//     using word_type = typename TypeParam::word_type;
//     auto const & fixture = this->fixture();
//     auto align_cfg = fixture.config;
//
//     std::vector database = fixture.sequence1;
//     std::vector query = fixture.sequence2;
//
//     auto alignment = edit_distance<word_type>(database, query, align_cfg);
//     auto trace_matrix = alignment.trace_matrix();
//     auto begin_coordinate = alignment.begin_coordinate();
//     auto end_coordinate = alignment.end_coordinate();
//
//     EXPECT_EQ(trace_matrix.cols(), database.size()+1);
//     EXPECT_EQ(trace_matrix.rows(), query.size()+1);
//     EXPECT_EQ(begin_coordinate.seq1_pos, fixture.begin_coordinate.seq1_pos);
//     EXPECT_EQ(begin_coordinate.seq2_pos, fixture.begin_coordinate.seq2_pos);
//     EXPECT_EQ(end_coordinate.seq1_pos, fixture.end_coordinate.seq1_pos);
//     EXPECT_EQ(end_coordinate.seq2_pos, fixture.end_coordinate.seq2_pos);
//     EXPECT_EQ(trace_matrix, fixture.trace_matrix);
//     EXPECT_EQ(alignment.score(), fixture.score);
//
//     auto && [gapped_database, gapped_query] = alignment.trace();
//     EXPECT_EQ(std::string{gapped_database | view::to_char}, fixture.gapped_sequence1);
//     EXPECT_EQ(std::string{gapped_query | view::to_char}, fixture.gapped_sequence2);
// }
//
// TYPED_TEST_P(global_affine_unbanded, trace)
// {
//     using word_type = typename TypeParam::word_type;
//     auto const & fixture = this->fixture();
//     auto align_cfg = fixture.config;
//
//     std::vector database = fixture.sequence1;
//     std::vector query = fixture.sequence2;
//
//     auto alignment = edit_distance<word_type>(database, query, align_cfg);
//
//     auto && [gapped_database, gapped_query] = alignment.trace();
//     EXPECT_EQ(std::string{gapped_database | view::to_char}, fixture.gapped_sequence1);
//     EXPECT_EQ(std::string{gapped_query | view::to_char}, fixture.gapped_sequence2);
// }

REGISTER_TYPED_TEST_CASE_P(global_affine_unbanded, score);

// work around a bug that you can't specify more than 50 template arguments to ::testing::types
INSTANTIATE_TYPED_TEST_CASE_P(global, global_affine_unbanded, global_affine_unbanded_types);
// INSTANTIATE_TYPED_TEST_CASE_P(semi_global, global_affine_unbanded, semi_global_edit_distance_unbanded_types);
// INSTANTIATE_TYPED_TEST_CASE_P(global_max_errors, global_affine_unbanded, global_edit_distance_max_errors_unbanded_types);
// INSTANTIATE_TYPED_TEST_CASE_P(semi_global_max_errors, global_affine_unbanded, semi_global_edit_distance_max_errors_unbanded_types);
