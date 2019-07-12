// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>

#include <seqan3/alignment/pairwise/execution/alignment_executor_two_way.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/range/view/view_all.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>
#include <seqan3/test/pretty_printing.hpp>

struct dummy_alignment
{

    template <typename first_seq_t, typename second_seq_t>
    constexpr auto operator()(size_t const, first_seq_t && first_seq, second_seq_t && second_seq) const
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

using namespace seqan3;

template <typename t>
struct alignment_executor_two_way_test : public ::testing::Test
{};

using testing_types = ::testing::Types<detail::execution_handler_sequential, detail::execution_handler_parallel>;
TYPED_TEST_CASE(alignment_executor_two_way_test, testing_types);

// Some globally defined test types
using seq_type = all_view<std::string &>;

inline static std::tuple single{std::string{"AACGTACGT"}, std::string{"ATCGTCCGT"}};
inline static std::vector<decltype(single)> collection{5, single};
inline static std::function<size_t(size_t const, seq_type, seq_type)> fn{dummy_alignment{}};

TYPED_TEST(alignment_executor_two_way_test, construction)
{
    using type = detail::alignment_executor_two_way<std::add_lvalue_reference_t<decltype(collection)>,
                                                    decltype(fn),
                                                    TypeParam>;

    EXPECT_FALSE(std::is_default_constructible_v<type>);
    EXPECT_FALSE(std::is_copy_constructible_v<type>);
    EXPECT_TRUE(std::is_move_constructible_v<type>);
    EXPECT_FALSE(std::is_copy_assignable_v<type>);
    EXPECT_TRUE(std::is_move_assignable_v<type>);
}

TYPED_TEST(alignment_executor_two_way_test, is_eof)
{
    using type = detail::alignment_executor_two_way<std::add_lvalue_reference_t<decltype(collection)>,
                                                    decltype(fn),
                                                    TypeParam>;
    type exec{collection, fn};
    EXPECT_FALSE(exec.is_eof());
}

TEST(alignment_executor_two_way, type_deduction)
{
    detail::alignment_executor_two_way exec{collection, fn};
    EXPECT_FALSE(exec.is_eof());
}

TYPED_TEST(alignment_executor_two_way_test, bump)
{
    using t = detail::alignment_executor_two_way<decltype(collection) &, decltype(fn), TypeParam>;
    t exec{collection, fn};

    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.bump()));
}

TYPED_TEST(alignment_executor_two_way_test, in_avail)
{
    using t = detail::alignment_executor_two_way<decltype(collection) &, decltype(fn), TypeParam>;
    t exec{collection, fn};

    EXPECT_EQ(exec.in_avail(), 0u);
    EXPECT_EQ(exec.bump().value(), 7u);
    if constexpr (std::Same<TypeParam, detail::execution_handler_parallel>)
        EXPECT_EQ(exec.in_avail(), 4u);
    else
        EXPECT_EQ(exec.in_avail(), 0u);
}

TYPED_TEST(alignment_executor_two_way_test, lvalue_single_view)
{
    auto v = ranges::view::single(single);
    using t = detail::alignment_executor_two_way<decltype(v), decltype(fn), TypeParam>;

    t exec{v, fn};
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.bump()));
}

TYPED_TEST(alignment_executor_two_way_test, rvalue_single_view)
{
    using t = detail::alignment_executor_two_way<decltype(ranges::view::single(single)), decltype(fn), TypeParam>;
    t exec{ranges::view::single(single), fn};
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.bump()));
}

TYPED_TEST(alignment_executor_two_way_test, lvalue_collection)
{
    using t = detail::alignment_executor_two_way<decltype(collection) &, decltype(fn), TypeParam>;
    t exec{collection, fn};
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.bump()));
}

TYPED_TEST(alignment_executor_two_way_test, rvalue_collection_view)
{
    using t = detail::alignment_executor_two_way<decltype(collection | view::persist), decltype(fn), TypeParam>;
    t exec{collection | view::persist, fn};
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_EQ(exec.bump().value(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.bump()));
}
