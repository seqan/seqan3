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
    constexpr void check_score(seqan3::detail::alignment_optimum<score_t> const &,
                               seqan3::detail::alignment_optimum<score_t> const &) const noexcept
    {}
};

template <typename cell_type>
class affine_gap_banded_policy_mock :
    public seqan3::detail::affine_gap_banded_policy<affine_gap_banded_policy_mock<cell_type>, cell_type>,
    public check_score_mock
{
public:
    using base_t = seqan3::detail::affine_gap_banded_policy<affine_gap_banded_policy_mock<cell_type>, cell_type>;

    using base_t::base_t;
    using base_t::compute_cell;
    using base_t::make_cache;
    using base_t::compute_first_band_cell;
};

TEST(affine_gap_banded_policy, construction)
{
    EXPECT_TRUE((std::is_default_constructible_v<affine_gap_banded_policy_mock<std::tuple<int, int, seqan3::detail::ignore_t>>>));
    EXPECT_TRUE((std::is_copy_constructible_v<affine_gap_banded_policy_mock<std::tuple<int, int, seqan3::detail::ignore_t>>>));
    EXPECT_TRUE((std::is_move_constructible_v<affine_gap_banded_policy_mock<std::tuple<int, int, seqan3::detail::ignore_t>>>));
    EXPECT_TRUE((std::is_copy_assignable_v<affine_gap_banded_policy_mock<std::tuple<int, int, seqan3::detail::ignore_t>>>));
    EXPECT_TRUE((std::is_move_assignable_v<affine_gap_banded_policy_mock<std::tuple<int, int, seqan3::detail::ignore_t>>>));
    EXPECT_TRUE((std::is_destructible_v<affine_gap_banded_policy_mock<std::tuple<int, int, seqan3::detail::ignore_t>>>));
}

TEST(affine_gap_banded_policy, make_cache)
{
    affine_gap_banded_policy_mock<std::tuple<int, int, seqan3::detail::ignore_t>> mock{};

    seqan3::gap_scheme scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}};
    auto cache = mock.make_cache(scheme);

    EXPECT_EQ(std::tuple_size_v<decltype(cache)>, 4u);

    using cell_t = std::tuple_element_t<0, decltype(cache)>;
    using gap_open_t = std::tuple_element_t<1, decltype(cache)>;
    using gap_extension_t = std::tuple_element_t<2, decltype(cache)>;
    using optimum_t = std::tuple_element_t<3, decltype(cache)>;

    EXPECT_TRUE((std::is_same_v<cell_t, std::tuple<int, int, seqan3::detail::ignore_t>>));
    EXPECT_TRUE((std::is_same_v<gap_open_t, int>));
    EXPECT_TRUE((std::is_same_v<gap_extension_t, int>));
    EXPECT_TRUE((std::is_same_v<optimum_t, seqan3::detail::alignment_optimum<int>>));

    EXPECT_EQ(std::get<1>(cache), -11);
    EXPECT_EQ(std::get<2>(cache), -1);
}

TEST(affine_gap_banded_policy, compute_first_band_cell)
{
    using namespace seqan3;
    using namespace seqan3::detail;

    affine_gap_banded_policy_mock<std::tuple<int, int, seqan3::detail::ignore_t>> mock{};

    gap_scheme scheme{gap_score{-1}, gap_open_score{-10}};
    auto cache = mock.make_cache(scheme);

    {  // max from diagonal
        std::tuple cell{std::tuple{0, -10, std::ignore}, std::tuple{-11, -20, std::ignore}};
        mock.compute_first_band_cell(std::make_tuple(std::ref(cell),
                                                     alignment_coordinate{column_index_type{3u}, row_index_type{5u}},
                                                     std::ignore),
                                     cache, 5);

        int first;
        int second;
        std::tie(first, second, std::ignore) = std::get<0>(cell);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{5, -10}));
        std::tie(first, second, std::ignore) = std::get<1>(cell);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{-11, -20}));
        std::tie(first, second, std::ignore) = std::get<0>(cache);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{0, -6}));
        EXPECT_EQ(std::get<1>(cache), -11);
        EXPECT_EQ(std::get<2>(cache), -1);
    }

    {  // max from horizontal
        std::tuple cell{std::tuple{0, -10, std::ignore}, std::tuple{-11, -20, std::ignore}};
        mock.compute_first_band_cell(std::make_tuple(std::ref(cell),
                                                     alignment_coordinate{column_index_type{3u}, row_index_type{5u}},
                                                     std::ignore),
                                     cache,
                                     -25);

        int first;
        int second;
        std::tie(first, second, std::ignore) = std::get<0>(cell);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{-20, -10}));
        std::tie(first, second, std::ignore) = std::get<1>(cell);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{-11, -20}));
        std::tie(first, second, std::ignore) = std::get<0>(cache);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{0, -31}));
        EXPECT_EQ(std::get<1>(cache), -11);
        EXPECT_EQ(std::get<2>(cache), -1);
    }

    {  // max from vertical - ignored in the computation.
        std::get<0>(cache) = std::tuple{0, 10, std::ignore};
        std::tuple cell{std::tuple{0, -10, std::ignore}, std::tuple{-11, -20, std::ignore}};
        mock.compute_first_band_cell(std::make_tuple(std::ref(cell),
                                                     alignment_coordinate{column_index_type{3u}, row_index_type{5u}},
                                                     std::ignore),
                                     cache,
                                     -10);

        int first;
        int second;
        std::tie(first, second, std::ignore) = std::get<0>(cell);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{-10, -10}));
        std::tie(first, second, std::ignore) = std::get<1>(cell);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{-11, -20}));
        std::tie(first, second, std::ignore) = std::get<0>(cache);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{0, -21}));
        EXPECT_EQ(std::get<1>(cache), -11);
        EXPECT_EQ(std::get<2>(cache), -1);
    }
}

TEST(affine_gap_banded_policy, compute_cell)
{
    using namespace seqan3;
    using namespace seqan3::detail;

    affine_gap_banded_policy_mock<std::tuple<int, int, seqan3::detail::ignore_t>> mock{};

    gap_scheme scheme{gap_score{-1}, gap_open_score{-10}};
    auto cache = mock.make_cache(scheme);

    { // max from diagonal
        std::get<0>(cache) = std::tuple{0, -4, std::ignore};
        std::tuple cell{std::tuple{0, -10, std::ignore}, std::tuple{-11, -20, std::ignore}};
        mock.compute_cell(std::make_tuple(std::ref(cell),
                                          alignment_coordinate{column_index_type{3u}, row_index_type{5u}},
                                          std::ignore),
                          cache,
                          5);
        int first;
        int second;
        std::tie(first, second, std::ignore) = std::get<0>(cell);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{5, -6}));
        std::tie(first, second, std::ignore) = std::get<1>(cell);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{-11, -20}));
        std::tie(first, second, std::ignore) = std::get<0>(cache);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{-6, -5}));
        EXPECT_EQ(std::get<1>(cache), -11);
        EXPECT_EQ(std::get<2>(cache), -1);
    }

    { // max from horizontal
        std::get<0>(cache) = std::tuple{0, -15, std::ignore};
        std::tuple cell{std::tuple{0, -10, std::ignore}, std::tuple{-11, -3, std::ignore}};
        mock.compute_cell(std::make_tuple(std::ref(cell),
                                          alignment_coordinate{column_index_type{3u}, row_index_type{5u}},
                                          std::ignore),
                          cache,
                          -10);

        int first;
        int second;
        std::tie(first, second, std::ignore) = std::get<0>(cell);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{-3, -4}));
        std::tie(first, second, std::ignore) = std::get<1>(cell);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{-11, -3}));
        std::tie(first, second, std::ignore) = std::get<0>(cache);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{-14, -14}));
        EXPECT_EQ(std::get<1>(cache), -11);
        EXPECT_EQ(std::get<2>(cache), -1);
    }

    {  // max from vertical
        std::get<0>(cache) = std::tuple{0, -3, std::ignore};
        std::tuple cell{std::tuple{0, -10, std::ignore}, std::tuple{-11, -4, std::ignore}};
        mock.compute_cell(std::make_tuple(std::ref(cell),
                                          alignment_coordinate{column_index_type{3u}, row_index_type{5u}},
                                          std::ignore),
                          cache,
                          -10);

        int first;
        int second;
        std::tie(first, second, std::ignore) = std::get<0>(cell);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{-3, -5}));
        std::tie(first, second, std::ignore) = std::get<1>(cell);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{-11, -4}));
        std::tie(first, second, std::ignore) = std::get<0>(cache);
        EXPECT_EQ((std::tie(first, second)), (std::tuple{-14, -4}));
        EXPECT_EQ(std::get<1>(cache), -11);
        EXPECT_EQ(std::get<2>(cache), -1);
    }
}
