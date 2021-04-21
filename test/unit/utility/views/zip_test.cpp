// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/std/algorithm>
#include <deque>
#include <iostream>
#include <list>
#include <seqan3/std/ranges>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/utility/views/zip.hpp>

// https://github.com/ericniebler/range-v3/issues/1480
TEST(zip_view, gcc10bug_rangev3_1480)
{
    // This regression test only checks if the respective code compiles.
    std::vector<char> const first_sequence{};
    std::vector<char> const second_sequence{};

    auto zip_view = seqan3::views::zip(first_sequence, second_sequence);
    std::ranges::for_each(zip_view, [&] ([[maybe_unused]] auto && value) {});
}
