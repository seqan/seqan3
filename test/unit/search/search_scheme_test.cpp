// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <type_traits>

#include <seqan3/search/detail/search_scheme_algorithm.hpp>

#include "helper_search_scheme.hpp"

#if SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY
using integral_t = uint16_t;
#else
using integral_t = uint8_t;
#endif // SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY

template <uint8_t min_error, uint8_t max_error, bool precomputed_scheme>
void error_distributions(auto & expected, auto & actual)
{
    if constexpr (precomputed_scheme)
    {
        auto const & oss{seqan3::detail::optimum_search_scheme<min_error, max_error>};
        seqan3::search_scheme_error_distribution(actual, oss);
        seqan3::search_scheme_error_distribution(
            expected,
            seqan3::trivial_search_scheme(min_error, max_error, oss.front().blocks()));
    }
    else
    {
        auto const & ss{seqan3::detail::compute_ss(min_error, max_error)};
        seqan3::search_scheme_error_distribution(actual, ss);
        seqan3::search_scheme_error_distribution(
            expected,
            seqan3::trivial_search_scheme(min_error, max_error, ss.front().blocks()));
    }
    std::sort(expected.begin(), expected.end());
    std::sort(actual.begin(), actual.end());
}

TEST(search_scheme_test, error_distribution_coverage_optimum_search_schemes)
{
    std::vector<std::vector<integral_t>> expected, actual;

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
    std::vector<std::vector<integral_t>> expected, actual;

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
    std::vector<std::vector<integral_t>> error_distributions;

    auto const & oss{seqan3::detail::optimum_search_scheme<min_error, max_error>};
    seqan3::search_scheme_error_distribution(error_distributions, oss);
    uint64_t size = error_distributions.size();
    std::sort(error_distributions.begin(), error_distributions.end());
    error_distributions.erase(std::unique(error_distributions.begin(), error_distributions.end()),
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
