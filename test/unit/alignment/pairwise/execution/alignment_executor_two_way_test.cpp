// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/single.hpp>
#include <range/v3/view/zip.hpp>

#include <seqan3/alignment/pairwise/execution/alignment_executor_two_way.hpp>
#include <seqan3/range/view/persist.hpp>

struct dummy_alignment
{

    template <typename first_seq_t, typename second_seq_t>
    constexpr auto operator()(first_seq_t && first_seq, second_seq_t && second_seq) const
    {
        size_t count = 0;
        ranges::for_each(ranges::view::zip(first_seq, second_seq), [&](auto && tpl)
        {
            auto && [v1, v2] = tpl;
            if (v1 == v2)
                ++count;
        });
        return count;
    }
};

// Some globally defined test types
inline static std::tuple single{std::string{"AACGTACGT"}, std::string{"ATCGTCCGT"}};
inline static std::vector<decltype(single)> collection{5, single};
inline static std::function<size_t(std::string const &, std::string const &)> fn{dummy_alignment{}};

using namespace seqan3;

TEST(alignment_executor_two_way, construction)
{
    using type = detail::alignment_executor_two_way<std::add_lvalue_reference_t<decltype(collection)>, decltype(fn)>;

    EXPECT_FALSE(std::is_default_constructible_v<type>);
    EXPECT_FALSE(std::is_copy_constructible_v<type>);
    EXPECT_TRUE(std::is_move_constructible_v<type>);
    EXPECT_FALSE(std::is_copy_assignable_v<type>);
    EXPECT_TRUE(std::is_move_assignable_v<type>);
}

TEST(alignment_executor_two_way, is_eof)
{
    using type = detail::alignment_executor_two_way<std::add_lvalue_reference_t<decltype(collection)>, decltype(fn)>;
    type exec{collection, fn};
    EXPECT_FALSE(exec.is_eof());
}

TEST(alignment_executor_two_way, type_deduction)
{
    detail::alignment_executor_two_way exec{collection, fn};
    EXPECT_FALSE(exec.is_eof());
}

TEST(alignment_executor_two_way, bump)
{
    detail::alignment_executor_two_way exec{collection, fn};
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.bump()));
}

TEST(alignment_executor_two_way, in_avail)
{
    detail::alignment_executor_two_way exec{collection, fn};
    EXPECT_EQ(exec.in_avail(), 0u);
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_EQ(exec.in_avail(), 0u);
}

TEST(alignment_executor_two_way, lvalue_single_view)
{
    auto v = ranges::view::single(single);
    detail::alignment_executor_two_way exec{v, fn};
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.bump()));
}

TEST(alignment_executor_two_way, rvalue_single_view)
{
    detail::alignment_executor_two_way exec{ranges::view::single(single), fn};
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.bump()));
}

TEST(alignment_executor_two_way, lvalue_collection)
{
    detail::alignment_executor_two_way exec{collection, fn};
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.bump()));
}

TEST(alignment_executor_two_way, rvalue_collection_view)
{
    detail::alignment_executor_two_way exec{collection | view::persist, fn};
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_EQ(exec.bump(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.bump()));
}
