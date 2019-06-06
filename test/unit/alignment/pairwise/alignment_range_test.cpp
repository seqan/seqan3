// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <optional>

#include <range/v3/view/iota.hpp>

#include <seqan3/alignment/pairwise/alignment_range.hpp>
#include <seqan3/range/view/single_pass_input.hpp>
#include <seqan3/range/shortcuts.hpp>

using namespace seqan3;

struct dummy_executor
{

    using value_type      = size_t;
    using reference       = value_type;
    using difference_type = std::ptrdiff_t;

    std::optional<size_t> bump()
    {
        auto it = begin(generator);
        if (it == end(generator))
        {
            return {};
        }
        else
        {
            std::optional<size_t> opt{*it};
            ++it;
            return opt;
        }
    }

private:

    detail::single_pass_input_view<decltype(ranges::view::iota(0u, 10u))> generator{ranges::view::iota(0u, 10u)};
};

TEST(alignment_range, concept_test)
{
    EXPECT_TRUE(std::ranges::InputRange<alignment_range<dummy_executor>>);
    EXPECT_FALSE(std::ranges::ForwardRange<alignment_range<dummy_executor>>);
}

TEST(alignment_range, construction)
{
    EXPECT_TRUE(std::is_default_constructible_v<alignment_range<dummy_executor>>);
    EXPECT_FALSE(std::is_copy_constructible_v<alignment_range<dummy_executor>>);
    EXPECT_TRUE(std::is_move_constructible_v<alignment_range<dummy_executor>>);
    EXPECT_FALSE(std::is_copy_assignable_v<alignment_range<dummy_executor>>);
    EXPECT_TRUE(std::is_move_assignable_v<alignment_range<dummy_executor>>);

    EXPECT_TRUE((std::is_constructible_v<alignment_range<dummy_executor>, dummy_executor>));
}

TEST(alignment_range, type_deduction)
{
    alignment_range rng{dummy_executor{}};
    EXPECT_TRUE((std::is_same_v<decltype(rng), alignment_range<dummy_executor>>));
}

TEST(alignment_range, begin)
{
    alignment_range rng{dummy_executor{}};
    auto it = rng.begin();
    EXPECT_EQ(*it, 0u);
}

TEST(alignment_range, end)
{
    alignment_range rng{dummy_executor{}};
    auto it = rng.end();
    EXPECT_FALSE(it == begin(rng));
    EXPECT_FALSE(begin(rng) == it);
}

TEST(alignment_range, iterable)
{
    alignment_range rng{dummy_executor{}};
    size_t sum = 0;
    for (size_t res : rng)
        sum += res;

    EXPECT_EQ(sum, 45u);
}

TEST(alignment_range, default_construction)
{
    alignment_range<dummy_executor> rng{};
    EXPECT_TRUE(begin(rng) == end(rng));
}
