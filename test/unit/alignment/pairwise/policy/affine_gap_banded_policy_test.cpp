// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/pairwise/policy/affine_gap_banded_policy.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>

struct check_score_mock
{
public:

    template <typename score_t>
    constexpr void check_score(score_t const &, score_t const &) const noexcept
    {}
};

template <typename cell_type>
class affine_gap_banded_policy_mock :
    public seqan3::detail::affine_gap_banded_policy<affine_gap_banded_policy_mock<cell_type>, cell_type>,
    public check_score_mock
{
public:
    using base_t = seqan3::detail::affine_gap_banded_policy<affine_gap_banded_policy_mock<cell_type>, cell_type>;

    using base_t::compute_cell;
    using base_t::make_cache;
    using base_t::compute_first_band_cell;
};

TEST(affine_gap_banded_policy, construction)
{
    EXPECT_TRUE((std::is_default_constructible_v<affine_gap_banded_policy_mock<std::tuple<int, int>>>));
    EXPECT_TRUE((std::is_copy_constructible_v<affine_gap_banded_policy_mock<std::tuple<int, int>>>));
    EXPECT_TRUE((std::is_move_constructible_v<affine_gap_banded_policy_mock<std::tuple<int, int>>>));
    EXPECT_TRUE((std::is_copy_assignable_v<affine_gap_banded_policy_mock<std::tuple<int, int>>>));
    EXPECT_TRUE((std::is_move_assignable_v<affine_gap_banded_policy_mock<std::tuple<int, int>>>));
    EXPECT_TRUE((std::is_destructible_v<affine_gap_banded_policy_mock<std::tuple<int, int>>>));
}

TEST(affine_gap_banded_policy, make_cache)
{
    affine_gap_banded_policy_mock<std::tuple<int, int>> mock{};

    seqan3::gap_scheme scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}};
    auto cache = mock.make_cache(scheme);

    EXPECT_EQ(std::tuple_size_v<decltype(cache)>, 4u);

    using cell_t = std::tuple_element_t<0, decltype(cache)>;
    using gap_open_t = std::tuple_element_t<1, decltype(cache)>;
    using gap_extension_t = std::tuple_element_t<2, decltype(cache)>;
    using optimum_t = std::tuple_element_t<3, decltype(cache)>;

    EXPECT_TRUE((std::is_same_v<cell_t, std::tuple<int, int>>));
    EXPECT_TRUE((std::is_same_v<gap_open_t, int>));
    EXPECT_TRUE((std::is_same_v<gap_extension_t, int>));
    EXPECT_TRUE((std::is_same_v<optimum_t, int>));

    EXPECT_EQ(std::get<1>(cache), -11);
    EXPECT_EQ(std::get<2>(cache), -1);
}

TEST(affine_gap_banded_policy, compute_first_band_cell)
{
    affine_gap_banded_policy_mock<std::tuple<int, int>> mock{};

    seqan3::gap_scheme scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}};
    auto cache = mock.make_cache(scheme);

    {  // max from diagonal
        std::tuple cell{std::tuple{0, -10}, std::tuple{-11, -20}};
        mock.compute_first_band_cell(cell, cache, 5);

        EXPECT_EQ(std::get<0>(cell), (std::tuple{5, -10}));
        EXPECT_EQ(std::get<1>(cell), (std::tuple{-11, -20}));
        EXPECT_EQ(std::get<0>(cache), (std::tuple{0, -6}));
        EXPECT_EQ(std::get<1>(cache), -11);
        EXPECT_EQ(std::get<2>(cache), -1);
    }

    {  // max from horizontal
        std::tuple cell{std::tuple{0, -10}, std::tuple{-11, -20}};
        mock.compute_first_band_cell(cell, cache, -25);

        EXPECT_EQ(std::get<0>(cell), (std::tuple{-20, -10}));
        EXPECT_EQ(std::get<1>(cell), (std::tuple{-11, -20}));
        EXPECT_EQ(std::get<0>(cache), (std::tuple{0, -31}));
        EXPECT_EQ(std::get<1>(cache), -11);
        EXPECT_EQ(std::get<2>(cache), -1);
    }

    {  // max from vertical - ignored in the computation.
        std::get<0>(cache) = std::tuple{0, 10};
        std::tuple cell{std::tuple{0, -10}, std::tuple{-11, -20}};
        mock.compute_first_band_cell(cell, cache, -10);

        EXPECT_EQ(std::get<0>(cell), (std::tuple{-10, -10}));
        EXPECT_EQ(std::get<1>(cell), (std::tuple{-11, -20}));
        EXPECT_EQ(std::get<0>(cache), (std::tuple{0, -21}));
        EXPECT_EQ(std::get<1>(cache), -11);
        EXPECT_EQ(std::get<2>(cache), -1);
    }
}

TEST(affine_gap_banded_policy, compute_cell)
{
    affine_gap_banded_policy_mock<std::tuple<int, int>> mock{};

    seqan3::gap_scheme scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}};
    auto cache = mock.make_cache(scheme);

    { // max from diagonal
        std::get<0>(cache) = std::tuple{0, -4};
        std::tuple cell{std::tuple{0, -10}, std::tuple{-11, -20}};
        mock.compute_cell(cell, cache, 5);

        EXPECT_EQ(std::get<0>(cell), (std::tuple{5, -6}));
        EXPECT_EQ(std::get<1>(cell), (std::tuple{-11, -20}));
        EXPECT_EQ(std::get<0>(cache), (std::tuple{-6, -5}));
        EXPECT_EQ(std::get<1>(cache), -11);
        EXPECT_EQ(std::get<2>(cache), -1);
    }

    { // max from horizontal
        std::get<0>(cache) = std::tuple{0, -15};
        std::tuple cell{std::tuple{0, -10}, std::tuple{-11, -3}};
        mock.compute_cell(cell, cache, -10);

        EXPECT_EQ(std::get<0>(cell), (std::tuple{-3, -4}));
        EXPECT_EQ(std::get<1>(cell), (std::tuple{-11, -3}));
        EXPECT_EQ(std::get<0>(cache), (std::tuple{-14, -14}));
        EXPECT_EQ(std::get<1>(cache), -11);
        EXPECT_EQ(std::get<2>(cache), -1);
    }

    {  // max from vertical
        std::get<0>(cache) = std::tuple{0, -3};
        std::tuple cell{std::tuple{0, -10}, std::tuple{-11, -4}};
        mock.compute_cell(cell, cache, -10);

        EXPECT_EQ(std::get<0>(cell), (std::tuple{-3, -5}));
        EXPECT_EQ(std::get<1>(cell), (std::tuple{-11, -4}));
        EXPECT_EQ(std::get<0>(cache), (std::tuple{-14, -4}));
        EXPECT_EQ(std::get<1>(cache), -11);
        EXPECT_EQ(std::get<2>(cache), -1);
    }
}
