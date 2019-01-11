// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <memory>

#include <range/v3/view/drop.hpp>
#include <range/v3/view/iota.hpp>

#include <seqan3/alignment/band/static_band.hpp>
#include <seqan3/alignment/pairwise/policy/banded_score_dp_matrix_policy.hpp>
#include <seqan3/range/view/take_exactly.hpp>

template <typename alloc_type>
class banded_score_dp_matrix_policy_mock :
    public seqan3::detail::banded_score_dp_matrix_policy<banded_score_dp_matrix_policy_mock<alloc_type>, alloc_type>
{
public:
    using base_t = seqan3::detail::banded_score_dp_matrix_policy<banded_score_dp_matrix_policy_mock<alloc_type>,
                                                                 alloc_type>;

    using base_t::base_t;

    // static variables
    using base_t::INF;

    // Member variables.
    using base_t::score_matrix;
    using base_t::dimension_first_range;
    using base_t::dimension_second_range;
    using base_t::current_column_index;
    using base_t::current_matrix_iter;
    using base_t::band_column_index;
    using base_t::band_row_index;

    // Member functions
    using base_t::second_range_begin_offset;
    using base_t::current_band_size;
    using base_t::next_column;
    using base_t::current_column;
    using base_t::allocate_matrix;
};

auto mock_factory(seqan3::static_band const & band = seqan3::static_band{seqan3::lower_bound{-3},
                                                                         seqan3::upper_bound{5}})
{
    banded_score_dp_matrix_policy_mock<std::allocator<std::tuple<int, int>>> mock{};

    std::string seq1 = "ACGTAGACTACTG";
    std::string seq2 = "ACGTAGACTACTGACGT";

    mock.allocate_matrix(seq1, seq2, band);
    return mock;
}

TEST(banded_score_dp_matrix_policy, construction)
{
    using mock_t = decltype(mock_factory());
    EXPECT_TRUE((std::is_default_constructible_v<mock_t>));
    EXPECT_TRUE((std::is_copy_constructible_v<mock_t>));
    EXPECT_TRUE((std::is_move_constructible_v<mock_t>));
    EXPECT_TRUE((std::is_copy_assignable_v<mock_t>));
    EXPECT_TRUE((std::is_move_assignable_v<mock_t>));
    EXPECT_TRUE((std::is_destructible_v<mock_t>));
}

TEST(banded_score_dp_matrix_policy, allocate_matrix)
{
    auto mock = mock_factory();

    EXPECT_EQ(mock.current_column_index, 0u);
    EXPECT_EQ(mock.score_matrix.size(), 10u);
    EXPECT_EQ(mock.dimension_first_range, 14u);
    EXPECT_EQ(mock.dimension_second_range, 18u);
    EXPECT_EQ(mock.band_column_index, 5u);
    EXPECT_EQ(mock.band_row_index, 3u);
    EXPECT_NE(mock.current_matrix_iter, seqan3::begin(mock.score_matrix));

    auto last_cell = *(--seqan3::end(mock.score_matrix));
    EXPECT_EQ(last_cell, (std::tuple{decltype(mock)::INF, decltype(mock)::INF}));
}

TEST(banded_score_dp_matrix_policy, next_column)
{
    auto mock = mock_factory();

    EXPECT_EQ(mock.current_column_index, 0u);
    EXPECT_NE(mock.current_matrix_iter, seqan3::begin(mock.score_matrix));

    mock.next_column();
    EXPECT_EQ(mock.current_column_index, 1u);
    EXPECT_NE(mock.current_matrix_iter, seqan3::begin(mock.score_matrix));
    mock.next_column();
    EXPECT_EQ(mock.current_column_index, 2u);
    EXPECT_NE(mock.current_matrix_iter, seqan3::begin(mock.score_matrix));
    mock.next_column();
    EXPECT_EQ(mock.current_column_index, 3u);
    EXPECT_NE(mock.current_matrix_iter, seqan3::begin(mock.score_matrix));
    mock.next_column();
    EXPECT_EQ(mock.current_column_index, 4u);
    EXPECT_NE(mock.current_matrix_iter, seqan3::begin(mock.score_matrix));
    mock.next_column();
    EXPECT_EQ(mock.current_column_index, 5u);
    EXPECT_EQ(mock.current_matrix_iter, seqan3::begin(mock.score_matrix));
    mock.next_column();
    EXPECT_EQ(mock.current_column_index, 6u);
    EXPECT_EQ(mock.current_matrix_iter, seqan3::begin(mock.score_matrix));
}

TEST(banded_score_dp_matrix_policy, current_band_size)
{
    using namespace seqan3;

    auto mock = mock_factory(static_band{lower_bound{-7}, upper_bound{5}});

    // Current band size is band_row_index + 1.
    for (auto s : ranges::view::iota(0u, 6u))
    {
        // Band increases by one as long as it is
        EXPECT_EQ(mock.current_band_size(), 8 + s);
        mock.next_column();
    }

    // After that the band size does not change until the end of second range is reached.
    for ([[maybe_unused]] auto s : ranges::view::iota(6u, 11u))
    {
        // Band increases by one as long as it is
        EXPECT_EQ(mock.current_band_size(), 13u);
        mock.next_column();
    }

    // When the band reaches the end it will be decreased by one
    for (auto s : ranges::view::iota(11u, 14u))
    {
        // Band increases by one as long as it is
        EXPECT_EQ(mock.current_band_size(), 13 - (s - 10u));
        mock.next_column();
    }
}

TEST(banded_score_dp_matrix_policy, current_column)
{
    using namespace seqan3;

    auto mock = mock_factory();

    auto col = mock.current_column() | detail::view_get_score_column;

    EXPECT_EQ(std::tuple_size_v<value_type_t<decltype(col)>>, 2u);
    EXPECT_EQ(seqan3::size(col), 4u);

    // Writing into the first value means writing into the second value
    for (auto && tpl : col)
    {
        std::get<0>(std::forward<decltype(tpl)>(tpl)) = std::tuple{-1, -1};
    }

    for (auto && tpl : col | view::take_exactly(3u))
        EXPECT_EQ(std::get<1>(std::forward<decltype(tpl)>(tpl)), (std::tuple{-1, -1}));

    EXPECT_EQ(std::get<1>(*std::ranges::prev(seqan3::end(col))),
              (std::tuple{decltype(mock)::INF, decltype(mock)::INF}));
}

TEST(banded_score_dp_matrix_policy, second_range_begin_offset)
{
    auto mock = mock_factory();

    // move to the first column behind the band position in first row.
    for (unsigned i = 0u; i < 6; ++i)
        mock.next_column();

    EXPECT_EQ(mock.second_range_begin_offset(), 0u);
    mock.next_column();
    EXPECT_EQ(mock.second_range_begin_offset(), 1u);
    mock.next_column();
    EXPECT_EQ(mock.second_range_begin_offset(), 2u);
    mock.next_column();
    EXPECT_EQ(mock.second_range_begin_offset(), 3u);
}

TEST(banded_score_dp_matrix_policy, band_touches_last_row)
{
    using namespace seqan3;

    auto mock = mock_factory(static_band{lower_bound{-7}, upper_bound{5}});

    for (unsigned i = 0; i < 10u; ++i)
    {
        EXPECT_FALSE(mock.band_touches_last_row());
        mock.next_column();
    }

    for (unsigned i = 10; i < 14u; ++i)
    {
        EXPECT_TRUE(mock.band_touches_last_row());
        mock.next_column();
    }
}

TEST(banded_score_dp_matrix_policy, trim_sequences)
{
    using namespace seqan3;

    banded_score_dp_matrix_policy_mock<std::allocator<std::tuple<int, int>>> mock{};
                      //0123456789
    std::string seq1 = "ACGTAGACTA";
    std::string seq2 = "ACGTAGACTA";

    {
        static_band band{lower_bound{-4}, upper_bound{4}};

        auto [t_seq1, t_seq2] = mock.trim_sequences(seq1, seq2, band);
        EXPECT_TRUE(ranges::equal(t_seq1, seq1));
        EXPECT_TRUE(ranges::equal(t_seq2, seq2));
        EXPECT_EQ(seqan3::size(t_seq1), seqan3::size(t_seq2));
    }

    {
        static_band band{lower_bound{3}, upper_bound{4}};

        auto [t_seq1, t_seq2] = mock.trim_sequences(seq1, seq2, band);
        EXPECT_TRUE(ranges::equal(t_seq1, seq1 | ranges::view::drop(2)));
        EXPECT_TRUE(ranges::equal(t_seq2, seq2 | view::take_exactly(7)));
    }

    {
        static_band band{lower_bound{-5}, upper_bound{-3}};

        auto [t_seq1, t_seq2] = mock.trim_sequences(seq1, seq2, band);
        EXPECT_TRUE(ranges::equal(t_seq1, seq1 | view::take_exactly(7)));
        EXPECT_TRUE(ranges::equal(t_seq2, seq2 | ranges::view::drop(2)));
    }
}
