// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <type_traits>

#include "helper_search_scheme.hpp"

#include <seqan3/search/detail/search_scheme_algorithm.hpp>

#include <gtest/gtest.h>

template <uint8_t min_error, uint8_t max_error, bool precomputed_scheme>
void error_distributions(auto & expected, auto & actual)
{
    if constexpr (precomputed_scheme)
    {
        auto const & oss{seqan3::detail::optimum_search_scheme<min_error, max_error>};
        seqan3::search_scheme_error_distribution(actual, oss);
        seqan3::search_scheme_error_distribution(expected, seqan3::trivial_search_scheme(min_error,
                                                                                         max_error,
                                                                                         oss.front().blocks()));
    }
    else
    {
        auto const & ss{seqan3::detail::compute_ss(min_error, max_error)};
        seqan3::search_scheme_error_distribution(actual, ss);
        seqan3::search_scheme_error_distribution(expected, seqan3::trivial_search_scheme(min_error,
                                                                                         max_error,
                                                                                         ss.front().blocks()));
    }
    std::sort(expected.begin(), expected.end());
    std::sort(actual.begin(), actual.end());
}

TEST(search_scheme_test, error_distribution_coverage_optimum_search_schemes)
{
    std::vector<std::vector<uint8_t> > expected, actual;

    error_distributions<0, 0, true>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<0, 1, true>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<1, 1, true>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<0, 2, true>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<1, 2, true>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<2, 2, true>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<0, 3, true>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<1, 3, true>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<2, 3, true>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<3, 3, true>(expected, actual);
    EXPECT_EQ(actual, expected);
}

TEST(search_scheme_test, error_distribution_coverage_computed_search_schemes)
{
    std::vector<std::vector<uint8_t> > expected, actual;

    error_distributions<0, 0, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<0, 1, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<1, 1, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<0, 2, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<1, 2, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<2, 2, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<0, 3, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<1, 3, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<2, 3, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<3, 3, false>(expected, actual);
    EXPECT_EQ(actual, expected);

    error_distributions<3, 5, false>(expected, actual);
    EXPECT_EQ(actual, expected);
    error_distributions<0, 6, false>(expected, actual);
    EXPECT_EQ(actual, expected);
    error_distributions<7, 7, false>(expected, actual);
    EXPECT_EQ(actual, expected);
}

template <uint8_t min_error, uint8_t max_error, bool precomputed_scheme>
bool check_disjoint_search_scheme()
{
    std::vector<std::vector<uint8_t> > error_distributions;

    auto const & oss{seqan3::detail::optimum_search_scheme<min_error, max_error>};
    seqan3::search_scheme_error_distribution(error_distributions, oss);
    uint64_t size = error_distributions.size();
    std::sort(error_distributions.begin(), error_distributions.end());
    error_distributions.erase(std::unique(error_distributions.begin(),
                                          error_distributions.end()),
                                          error_distributions.end());
    return size == error_distributions.size();
}

TEST(search_scheme_test, error_distribution_disjoint_optimum_search_schemes)
{
    bool ret;

    ret = check_disjoint_search_scheme<0, 0, false>();
    EXPECT_TRUE(ret);
    ret = check_disjoint_search_scheme<0, 1, false>();
    EXPECT_TRUE(ret);
    ret = check_disjoint_search_scheme<0, 2, false>();
    EXPECT_TRUE(ret);
    ret = check_disjoint_search_scheme<0, 3, false>();
    EXPECT_TRUE(ret);
}
