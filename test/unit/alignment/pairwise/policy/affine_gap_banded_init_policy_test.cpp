// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/band/static_band.hpp>
#include <seqan3/alignment/pairwise/policy/affine_gap_banded_init_policy.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>

class affine_gap_banded_init_policy_mock :
    public seqan3::detail::affine_gap_banded_init_policy<affine_gap_banded_init_policy_mock>
{
public:
    using base_t = seqan3::detail::affine_gap_banded_init_policy<affine_gap_banded_init_policy_mock>;

    using base_t::init_origin_cell;
    using base_t::init_column_cell;
    using base_t::init_row_cell;
    using base_t::balance_leading_gaps;
};

TEST(affine_gap_banded_init_policy, construction)
{
    EXPECT_TRUE(std::is_default_constructible_v<affine_gap_banded_init_policy_mock>);
    EXPECT_TRUE(std::is_copy_constructible_v<affine_gap_banded_init_policy_mock>);
    EXPECT_TRUE(std::is_move_constructible_v<affine_gap_banded_init_policy_mock>);
    EXPECT_TRUE(std::is_copy_assignable_v<affine_gap_banded_init_policy_mock>);
    EXPECT_TRUE(std::is_move_assignable_v<affine_gap_banded_init_policy_mock>);
    EXPECT_TRUE(std::is_destructible_v<affine_gap_banded_init_policy_mock>);
}

TEST(affine_gap_banded_init_policy, init_origin_cell)
{
    std::tuple cell{std::tuple{0, 0}, std::tuple{0, 0}};
    std::tuple cache{std::tuple{0, 0}, -10, -1};

    affine_gap_banded_init_policy_mock mock{};

    mock.init_origin_cell(cell, cache);

    EXPECT_EQ(std::get<0>(cell), (std::tuple{0, -10}));
    EXPECT_EQ(std::get<1>(cell), (std::tuple{0, 0}));
    EXPECT_EQ(std::get<0>(cache), (std::tuple{0, -10}));
    EXPECT_EQ(std::get<1>(cache), -10);
    EXPECT_EQ(std::get<2>(cache), -1);
}

TEST(affine_gap_banded_init_policy, init_column_cell)
{
    std::tuple cell{std::tuple{0, -10}, std::tuple{0, 0}};
    std::tuple cache{std::tuple{0, -10}, -10, -1};

    affine_gap_banded_init_policy_mock mock{};

    mock.init_column_cell(cell, cache);

    EXPECT_EQ(std::get<0>(cell), (std::tuple{-10, -20}));
    EXPECT_EQ(std::get<1>(cell), (std::tuple{0, 0}));
    EXPECT_EQ(std::get<0>(cache), (std::tuple{0, -11}));
    EXPECT_EQ(std::get<1>(cache), -10);
    EXPECT_EQ(std::get<2>(cache), -1);
}

TEST(affine_gap_banded_init_policy, init_row_cell)
{
    std::tuple cell{std::tuple{0, 0}, std::tuple{0, -10}};
    std::tuple cache{std::tuple{0, 0}, -10, -1};

    affine_gap_banded_init_policy_mock mock{};

    mock.init_row_cell(cell, cache);

    EXPECT_EQ(std::get<0>(cell), (std::tuple{-10, -11}));
    EXPECT_EQ(std::get<1>(cell), (std::tuple{0, -10}));
    EXPECT_EQ(std::get<0>(cache), (std::tuple{0, -20}));
    EXPECT_EQ(std::get<1>(cache), -10);
    EXPECT_EQ(std::get<2>(cache), -1);
}

// TODO Templatize
TEST(affine_gap_banded_init_policy, balance_leading_gaps)
{
    using namespace seqan3;
    //
    static_band band{lower_bound{-3}, upper_bound{3}};
    gap_scheme scheme{gap_score{-1}, gap_open_score{-10}};
    affine_gap_banded_init_policy_mock mock{};
    //
    int total = 0;

    mock.balance_leading_gaps(total, band, scheme);
    EXPECT_EQ(total, 0);

    band.lower_bound = -4;
    band.upper_bound = -3;

    mock.balance_leading_gaps(total, band, scheme);
    EXPECT_EQ(total, -13);

    band.lower_bound = 4;
    band.upper_bound = 10;

    mock.balance_leading_gaps(total, band, scheme);
    EXPECT_EQ(total, -27);
}
