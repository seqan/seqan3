// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides test utilities for std::ranges::range types.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>

#include <seqan3/core/debug_stream/detail/to_string.hpp>
#include <seqan3/core/debug_stream/range.hpp>
#include <seqan3/test/pretty_printing.hpp>

namespace seqan3::test
{
#define EXPECT_RANGE_EQ(val1, val2) EXPECT_PRED_FORMAT2(::seqan3::test::expect_range_eq{}, val1, val2);

struct expect_range_eq
{
    template <typename rng_t>
    auto copy_range(rng_t && rng)
    {
        using value_t = std::ranges::range_value_t<rng_t>;
        std::vector<value_t> rng_copy{};
        std::ranges::copy(std::forward<rng_t>(rng), std::back_inserter(rng_copy));
        return rng_copy;
    }

    template <std::ranges::range lhs_t, std::ranges::range rhs_t>
    ::testing::AssertionResult
    operator()(char const * lhs_expression, char const * rhs_expression, lhs_t && lhs, rhs_t && rhs)
    {
        std::vector lhs_copy = copy_range(std::forward<lhs_t>(lhs));
        std::vector rhs_copy = copy_range(std::forward<rhs_t>(rhs));

        if (std::ranges::equal(lhs_copy, rhs_copy))
            return ::testing::AssertionSuccess();

        return ::testing::internal::CmpHelperEQFailure(lhs_expression, rhs_expression, lhs_copy, rhs_copy);
    }
};

} // namespace seqan3::test
