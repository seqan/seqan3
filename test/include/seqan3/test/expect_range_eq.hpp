// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides test utilities for std::ranges::range types.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <gtest/gtest.h>

#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

#include <seqan3/core/detail/debug_stream_range.hpp>
#include <seqan3/core/detail/to_string.hpp>
#include <seqan3/test/pretty_printing.hpp>

namespace seqan3::test
{
#define EXPECT_RANGE_EQ(val1, val2) \
    EXPECT_PRED_FORMAT2(::seqan3::test::expect_range_eq{}, val1, val2);

struct expect_range_eq
{
    template <typename rng_t>
    auto copy_range(rng_t && rng)
    {
        using value_t = std::ranges::range_value_t<rng_t>;
        std::vector<value_t> rng_copy{};
        std::ranges::copy(rng, std::ranges::back_inserter(rng_copy));
        return rng_copy;
    }

    template <std::ranges::range lhs_t, std::ranges::range rhs_t>
    ::testing::AssertionResult operator()(char const * lhs_expression, char const * rhs_expression,
                                          lhs_t && lhs, rhs_t && rhs)
    {
        std::vector lhs_copy = copy_range(lhs);
        std::vector rhs_copy = copy_range(rhs);

        if (std::ranges::equal(lhs_copy, rhs_copy))
            return ::testing::AssertionSuccess();

        return ::testing::internal::CmpHelperEQFailure(lhs_expression, rhs_expression, lhs_copy, rhs_copy);
    }
};

} // namespace seqan3::test
