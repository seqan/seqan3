// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>
#include <memory>

#include <gtest/gtest.h>

#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/range/views/detail.hpp>

// The general capabilities of adaptor_base and derivates are tested thoroughly by the different views
// this file checks the correct memory behaviour in regard to storing the elements
// (hold and pass references if possible; for values move-in/out if possible)

struct copy_counter
{
    size_t copy_count = 0;
    size_t move_count = 0;

    copy_counter() = default;
    copy_counter(copy_counter const & rhs)
    {
        copy_count = rhs.copy_count + 1;
        move_count = rhs.move_count;
    }

    copy_counter(copy_counter && rhs)
    {
        copy_count = rhs.copy_count;
        move_count = rhs.move_count + 1;
    }

    copy_counter & operator=(copy_counter const &) =  delete;
    copy_counter & operator=(copy_counter &&) =  delete;
};

struct adaptor_base_type_checker :
    seqan3::detail::adaptor_base<adaptor_base_type_checker,
                                 copy_counter, copy_counter const, copy_counter &, copy_counter const &>
{
    using base_t = seqan3::detail::adaptor_base<adaptor_base_type_checker,
                                                copy_counter, copy_counter const, copy_counter &, copy_counter const &>;
    using base_t::base_t;

    template <typename urng_t, typename one_t, typename two_t, typename three_t, typename four_t>
    static std::tuple<one_t, two_t, three_t, four_t>
    impl(urng_t &&, one_t && one, two_t && two, three_t && three, four_t && four)
    {
        return { std::forward<one_t>(one),
                 std::forward<two_t>(two),
                 std::forward<three_t>(three),
                 std::forward<four_t>(four) };
    }
};

TEST(arg_ownership, lval_adaptor)
{
    copy_counter c3, c4;

    adaptor_base_type_checker a{copy_counter{}, copy_counter{}, c3, c4};

    std::vector<int> vec;

    auto f = vec | a;

    EXPECT_TRUE((std::same_as<decltype(f),
                              std::tuple<copy_counter, copy_counter const, copy_counter &, copy_counter const &>>));

    // In general three operations happen:
    // 1. out of constructor, into storage tuple
    // 2. out of storage tuple, into impl()
    // 3. from impl(), into return tuple
    EXPECT_EQ(std::get<0>(f).copy_count, 1ul);  // 2. because needs to stay
    EXPECT_EQ(std::get<0>(f).move_count, 2ul);  // 1. and 3.

    EXPECT_EQ(std::get<1>(f).copy_count, 3ul);
    EXPECT_EQ(std::get<1>(f).move_count, 0ul);

    EXPECT_EQ(std::get<2>(f).copy_count, 0ul);
    EXPECT_EQ(std::get<2>(f).move_count, 0ul);
    EXPECT_EQ(c3.copy_count, 0ul);
    EXPECT_EQ(c3.move_count, 0ul);

    EXPECT_EQ(std::get<3>(f).copy_count, 0ul);
    EXPECT_EQ(std::get<3>(f).move_count, 0ul);
    EXPECT_EQ(c4.copy_count, 0ul);
    EXPECT_EQ(c4.move_count, 0ul);
}

TEST(arg_ownership, const_lval_adaptor)
{
    copy_counter c3, c4;

    adaptor_base_type_checker const a{copy_counter{}, copy_counter{}, c3, c4};

    std::vector<int> vec;

    auto f = vec | a;

    EXPECT_TRUE((std::same_as<decltype(f),
                              std::tuple<copy_counter, copy_counter const, copy_counter &, copy_counter const &>>));

    EXPECT_EQ(std::get<0>(f).copy_count, 1ul);
    EXPECT_EQ(std::get<0>(f).move_count, 2ul);

    EXPECT_EQ(std::get<1>(f).copy_count, 3ul);
    EXPECT_EQ(std::get<1>(f).move_count, 0ul);

    EXPECT_EQ(std::get<2>(f).copy_count, 0ul);
    EXPECT_EQ(std::get<2>(f).move_count, 0ul);
    EXPECT_EQ(c3.copy_count, 0ul);
    EXPECT_EQ(c3.move_count, 0ul);

    EXPECT_EQ(std::get<3>(f).copy_count, 0ul);
    EXPECT_EQ(std::get<3>(f).move_count, 0ul);
    EXPECT_EQ(c4.copy_count, 0ul);
    EXPECT_EQ(c4.move_count, 0ul);
}

TEST(arg_ownership, rval_adaptor)
{
    copy_counter c3, c4;

    adaptor_base_type_checker a{copy_counter{}, copy_counter{}, c3, c4};

    std::vector<int> vec;

    auto f = vec | std::move(a);

    EXPECT_TRUE((std::same_as<decltype(f),
                              std::tuple<copy_counter, copy_counter const, copy_counter &, copy_counter const &>>));

    EXPECT_EQ(std::get<0>(f).copy_count, 0ul); // moved out of storage, too, because temporary
    EXPECT_EQ(std::get<0>(f).move_count, 3ul);

    EXPECT_EQ(std::get<1>(f).copy_count, 3ul);
    EXPECT_EQ(std::get<1>(f).move_count, 0ul);

    EXPECT_EQ(std::get<2>(f).copy_count, 0ul);
    EXPECT_EQ(std::get<2>(f).move_count, 0ul);
    EXPECT_EQ(c3.copy_count, 0ul);
    EXPECT_EQ(c3.move_count, 0ul);

    EXPECT_EQ(std::get<3>(f).copy_count, 0ul);
    EXPECT_EQ(std::get<3>(f).move_count, 0ul);
    EXPECT_EQ(c4.copy_count, 0ul);
    EXPECT_EQ(c4.move_count, 0ul);
}

TEST(arg_ownership, const_rval_adaptor)
{
    copy_counter c3, c4;

    adaptor_base_type_checker const a{copy_counter{}, copy_counter{}, c3, c4};

    std::vector<int> vec;

    auto f = vec | std::move(a);

    EXPECT_TRUE((std::same_as<decltype(f),
                              std::tuple<copy_counter, copy_counter const, copy_counter &, copy_counter const &>>));

    EXPECT_EQ(std::get<0>(f).copy_count, 1ul);
    EXPECT_EQ(std::get<0>(f).move_count, 2ul);

    EXPECT_EQ(std::get<1>(f).copy_count, 3ul);
    EXPECT_EQ(std::get<1>(f).move_count, 0ul);

    EXPECT_EQ(std::get<2>(f).copy_count, 0ul);
    EXPECT_EQ(std::get<2>(f).move_count, 0ul);
    EXPECT_EQ(c3.copy_count, 0ul);
    EXPECT_EQ(c3.move_count, 0ul);

    EXPECT_EQ(std::get<3>(f).copy_count, 0ul);
    EXPECT_EQ(std::get<3>(f).move_count, 0ul);
    EXPECT_EQ(c4.copy_count, 0ul);
    EXPECT_EQ(c4.move_count, 0ul);
}

template <typename t>
struct dummy_view
{};

TEST(adaptor_combination, constexpr_combine)
{
    constexpr auto adaptor1 = seqan3::detail::adaptor_for_view_without_args<dummy_view>{};
    constexpr auto adaptor2 = seqan3::detail::adaptor_for_view_without_args<dummy_view>{};

    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(adaptor1 | adaptor2)));
}
