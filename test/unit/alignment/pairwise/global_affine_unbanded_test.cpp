// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>


#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>

#include <seqan3/range/view/to_char.hpp>

#include "fixture/global_affine_unbanded.hpp"
#include "fixture/semi_global_affine_unbanded.hpp"

using namespace seqan3;
using namespace seqan3::detail;
using namespace seqan3::test::alignment::fixture;

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
        param<&global::affine::unbanded::dna4_02>
    >;

using semi_global_affine_unbanded_types
    = ::testing::Types<
        param<&semi_global::affine::unbanded::dna4_01>,
        param<&semi_global::affine::unbanded::dna4_02>,
        param<&semi_global::affine::unbanded::dna4_03>,
        param<&semi_global::affine::unbanded::dna4_04>
    >;

TYPED_TEST_P(global_affine_unbanded, score)
{
    auto const & fixture = this->fixture();
    auto align_cfg = fixture.config | align_cfg::result{with_score};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = align_pairwise(std::tie(database, query), align_cfg);

    EXPECT_EQ((*std::ranges::begin(alignment)).score(), fixture.score);
}

TYPED_TEST_P(global_affine_unbanded, end_position)
{
    auto const & fixture = this->fixture();
    auto align_cfg = fixture.config | align_cfg::result{with_back_coordinate};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = align_pairwise(std::tie(database, query), align_cfg);

    auto res = *std::ranges::begin(alignment);
    EXPECT_EQ(res.score(), fixture.score);
    EXPECT_EQ(res.back_coordinate(), fixture.back_coordinate);
}

TYPED_TEST_P(global_affine_unbanded, begin_position)
{
    auto const & fixture = this->fixture();
    auto align_cfg = fixture.config | align_cfg::result{with_front_coordinate};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = align_pairwise(std::tie(database, query), align_cfg);

    auto res = *std::ranges::begin(alignment);
    EXPECT_EQ(res.score(), fixture.score);
    EXPECT_EQ(res.front_coordinate(), fixture.front_coordinate);
    EXPECT_EQ(res.back_coordinate(), fixture.back_coordinate);
}

TYPED_TEST_P(global_affine_unbanded, trace)
{
    auto const & fixture = this->fixture();
    auto align_cfg = fixture.config | align_cfg::result{with_alignment};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = align_pairwise(std::tie(database, query), align_cfg);

    auto res = *std::ranges::begin(alignment);
    EXPECT_EQ(res.score(), fixture.score);
    EXPECT_EQ(res.front_coordinate(), fixture.front_coordinate);
    EXPECT_EQ(res.back_coordinate(), fixture.back_coordinate);
    EXPECT_TRUE(ranges::equal(get<0>(res.alignment()) | view::to_char, fixture.aligned_sequence1));
    EXPECT_TRUE(ranges::equal(get<1>(res.alignment()) | view::to_char, fixture.aligned_sequence2));
}

REGISTER_TYPED_TEST_CASE_P(global_affine_unbanded, score, end_position, begin_position, trace);

INSTANTIATE_TYPED_TEST_CASE_P(global, global_affine_unbanded, global_affine_unbanded_types);
INSTANTIATE_TYPED_TEST_CASE_P(semi_global, global_affine_unbanded, semi_global_affine_unbanded_types);
