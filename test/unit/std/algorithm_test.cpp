// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/algorithm>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>
#include <vector>

#include <seqan3/core/platform.hpp>

// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95578
TEST(algorithm_test, gcc95578)
{
    std::vector<int> v{};
    auto && rng = v | std::views::take_while([](auto &&){return true;});

    std::vector<int> rng_copy{};
    std::ranges::copy(rng, std::cpp20::back_inserter(rng_copy));
}
