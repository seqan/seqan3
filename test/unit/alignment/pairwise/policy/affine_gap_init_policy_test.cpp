// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/pairwise/detail/alignment_algorithm_state.hpp>
#include <seqan3/alignment/matrix/detail/alignment_score_matrix_one_column_banded.hpp>
#include <seqan3/alignment/matrix/detail/alignment_score_matrix_one_column.hpp>
#include <seqan3/alignment/matrix/detail/alignment_trace_matrix_full_banded.hpp>
#include <seqan3/alignment/matrix/detail/alignment_trace_matrix_full.hpp>
#include <seqan3/alignment/pairwise/policy/affine_gap_init_policy.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/std/ranges>

using seqan3::operator""_dna4;

class affine_gap_init_policy_mock :
    public seqan3::detail::affine_gap_init_policy<affine_gap_init_policy_mock>
{
public:
    using base_t = seqan3::detail::affine_gap_init_policy<affine_gap_init_policy_mock>;

    using base_t::base_t;
    using base_t::init_origin_cell;
    using base_t::init_column_cell;
    using base_t::init_row_cell;

    template <typename cell_t, typename state_t>
    void check_score_of_cell(cell_t, state_t) const
    {
        // no-op
    }
};

template <typename test_types>
class affine_gap_init_fixture : public ::testing::Test
{
public:
    using score_matrix_t = std::tuple_element_t<0, test_types>;
    using trace_matrix_t = std::tuple_element_t<1, test_types>;

    using score_matrix_iter_t = typename score_matrix_t::iterator;
    using trace_matrix_iter_t = typename trace_matrix_t::iterator;

    static constexpr bool with_trace =
        !seqan3::detail::decays_to_ignore_v<
            std::remove_reference_t<decltype(std::declval<typename trace_matrix_t::value_type>().current)>>;

    void SetUp()
    {
        if constexpr (std::tuple_element_t<2, test_types>::value)
        {
            score_matrix = score_matrix_t{"ACGT"_dna4, "ACGT"_dna4, seqan3::static_band{seqan3::lower_bound{-2},
                                                                                        seqan3::upper_bound{2}}};
            trace_matrix = trace_matrix_t{"ACGT"_dna4, "ACGT"_dna4, seqan3::static_band{seqan3::lower_bound{-2},
                                                                                        seqan3::upper_bound{2}}};
        }
        else
        {
            score_matrix = score_matrix_t{"ACGT"_dna4, "ACGT"_dna4};
            trace_matrix = trace_matrix_t{"ACGT"_dna4, "ACGT"_dna4};
        }
        score_matrix_iter = score_matrix.begin();
        trace_matrix_iter = trace_matrix.begin();
        state.gap_open_score = -10;
        state.gap_extension_score = -1;
    }

    auto column()
    {
        return seqan3::views::zip(*score_matrix_iter, *trace_matrix_iter);
    }

    affine_gap_init_policy_mock mock{};
    seqan3::detail::alignment_algorithm_state<int> state{};

    score_matrix_t score_matrix{};
    trace_matrix_t trace_matrix{};

    score_matrix_iter_t score_matrix_iter{};
    trace_matrix_iter_t trace_matrix_iter{};
};

using testing_types = ::testing::Types<
    std::tuple<seqan3::detail::alignment_score_matrix_one_column<int>,
               seqan3::detail::alignment_trace_matrix_full<seqan3::detail::trace_directions>,
               std::false_type>,
    std::tuple<seqan3::detail::alignment_score_matrix_one_column<int>,
               seqan3::detail::alignment_trace_matrix_full<seqan3::detail::trace_directions, true>,
               std::false_type>,
    std::tuple<seqan3::detail::alignment_score_matrix_one_column_banded<int>,
               seqan3::detail::alignment_trace_matrix_full_banded<seqan3::detail::trace_directions>,
               std::true_type>,
    std::tuple<seqan3::detail::alignment_score_matrix_one_column_banded<int>,
               seqan3::detail::alignment_trace_matrix_full_banded<seqan3::detail::trace_directions, true>,
               std::true_type>
    >;

TYPED_TEST_SUITE(affine_gap_init_fixture, testing_types, );

TYPED_TEST(affine_gap_init_fixture, init_origin_cell)
{
    auto zip_column = this->column();
    auto it = zip_column.begin();

    this->mock.init_origin_cell(*it, this->state);

    auto [score_cell, trace_cell] = *it;

    EXPECT_EQ(score_cell.current, 0);
    EXPECT_EQ(score_cell.up, -10);
    EXPECT_EQ(score_cell.w_left, -10);

    if constexpr (TestFixture::with_trace)
    {
        EXPECT_EQ(trace_cell.current, seqan3::detail::trace_directions::none);
        EXPECT_EQ(trace_cell.up, seqan3::detail::trace_directions::up_open);
        EXPECT_EQ(trace_cell.w_left, seqan3::detail::trace_directions::left_open);
    }
}

TYPED_TEST(affine_gap_init_fixture, init_column_cell)
{
    auto zip_column = this->column();
    auto it = zip_column.begin();

    this->mock.init_origin_cell(*it, this->state);
    ++it;
    this->mock.init_column_cell(*it, this->state);

    auto [score_cell, trace_cell] = *it;

    EXPECT_EQ(score_cell.current, -10);
    EXPECT_EQ(score_cell.up, -11);
    EXPECT_EQ(score_cell.w_left, -20);

    if constexpr (TestFixture::with_trace)
    {
        EXPECT_EQ(trace_cell.current, seqan3::detail::trace_directions::up_open);
        EXPECT_EQ(trace_cell.up, seqan3::detail::trace_directions::up);
        EXPECT_EQ(trace_cell.w_left, seqan3::detail::trace_directions::left_open);
    }
}

TYPED_TEST(affine_gap_init_fixture, init_row_cell)
{
    auto zip_column = this->column();
    auto it = zip_column.begin();

    this->mock.init_origin_cell(*it, this->state);
    ++this->score_matrix_iter;
    ++this->trace_matrix_iter;

    zip_column = this->column();
    it = zip_column.begin();
    this->mock.init_row_cell(*it, this->state);
    auto [score_cell, trace_cell] = *it;

    EXPECT_EQ(score_cell.current, -10);
    EXPECT_EQ(score_cell.up, -20);
    EXPECT_EQ(score_cell.w_left, -11);

    if constexpr (TestFixture::with_trace)
    {
        EXPECT_EQ(trace_cell.current, seqan3::detail::trace_directions::left_open);
        EXPECT_EQ(trace_cell.up, seqan3::detail::trace_directions::up_open);
        EXPECT_EQ(trace_cell.w_left, seqan3::detail::trace_directions::left);
    }
}
