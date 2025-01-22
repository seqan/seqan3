// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iostream>
#include <memory>
#include <ranges>

#include <seqan3/core/range/detail/adaptor_base.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

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

    copy_counter & operator=(copy_counter const &) = delete;
    copy_counter & operator=(copy_counter &&) = delete;
};

struct adaptor_base_type_checker :
    seqan3::detail::
        adaptor_base<adaptor_base_type_checker, copy_counter, copy_counter const, copy_counter &, copy_counter const &>
{
    using base_t = seqan3::detail::
        adaptor_base<adaptor_base_type_checker, copy_counter, copy_counter const, copy_counter &, copy_counter const &>;
    using base_t::base_t;

    template <typename urng_t, typename one_t, typename two_t, typename three_t, typename four_t>
    static std::tuple<one_t, two_t, three_t, four_t>
    impl(urng_t &&, one_t && one, two_t && two, three_t && three, four_t && four)
    {
        return {std::forward<one_t>(one),
                std::forward<two_t>(two),
                std::forward<three_t>(three),
                std::forward<four_t>(four)};
    }
};

TEST(arg_ownership, lval_adaptor)
{
    copy_counter c3, c4;

    adaptor_base_type_checker a{copy_counter{}, copy_counter{}, c3, c4};

    std::vector<int> vec;

    auto f = vec | a;

    EXPECT_TRUE((
        std::same_as<decltype(f), std::tuple<copy_counter, copy_counter const, copy_counter &, copy_counter const &>>));

    // In general three operations happen:
    // 1. out of constructor, into storage tuple
    // 2. out of storage tuple, into impl()
    // 3. from impl(), into return tuple
    EXPECT_EQ(std::get<0>(f).copy_count, 1ul); // 2. because needs to stay
    EXPECT_EQ(std::get<0>(f).move_count, 2ul); // 1. and 3.

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

    EXPECT_TRUE((
        std::same_as<decltype(f), std::tuple<copy_counter, copy_counter const, copy_counter &, copy_counter const &>>));

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

    EXPECT_TRUE((
        std::same_as<decltype(f), std::tuple<copy_counter, copy_counter const, copy_counter &, copy_counter const &>>));

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

    EXPECT_TRUE((
        std::same_as<decltype(f), std::tuple<copy_counter, copy_counter const, copy_counter &, copy_counter const &>>));

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

struct take_adaptor_closure : public seqan3::detail::adaptor_base<take_adaptor_closure, size_t>
{
    using base_t = seqan3::detail::adaptor_base<take_adaptor_closure, size_t>;
    using base_t::base_t;

    template <typename urng_t>
    static auto impl(urng_t && urng, size_t size)
    {
        using difference_t = std::ranges::range_difference_t<urng_t>;
        return std::views::take(std::forward<urng_t>(urng), static_cast<difference_t>(size));
    }
};

struct drop_adaptor_closure : public seqan3::detail::adaptor_base<drop_adaptor_closure, size_t>
{
    using base_t = seqan3::detail::adaptor_base<drop_adaptor_closure, size_t>;
    using base_t::base_t;

    template <typename urng_t>
    static auto impl(urng_t && urng, size_t size)
    {
        using difference_t = std::ranges::range_difference_t<urng_t>;
        return std::views::drop(std::forward<urng_t>(urng), static_cast<difference_t>(size));
    }
};

TEST(adaptor_base, function_object)
{
    constexpr take_adaptor_closure take1{1};
    constexpr drop_adaptor_closure drop1{1};
    std::vector<int> vec{0, 1, 2, 3, 4, 5};

    EXPECT_SAME_TYPE(decltype(take1(vec)), decltype(std::views::take(vec, 1)));
    EXPECT_RANGE_EQ(take1(vec), (std::vector<int>{0}));

    EXPECT_SAME_TYPE(decltype(drop1(vec)), decltype(std::views::drop(vec, 1)));
    EXPECT_RANGE_EQ(drop1(vec), (std::vector<int>{1, 2, 3, 4, 5}));
}

TEST(adaptor_base, pipe_range)
{
    constexpr take_adaptor_closure take1{1};
    std::vector<int> vec{0, 1, 2, 3, 4, 5};

    EXPECT_SAME_TYPE(decltype(vec | take1), decltype(std::views::take(vec, 1)));
    EXPECT_RANGE_EQ(vec | take1, (std::vector<int>{0}));
}

TEST(adaptor_base, pipe_same_adaptor)
{
    constexpr take_adaptor_closure take1{1};
    constexpr take_adaptor_closure take3{3};
    std::vector<int> vec{0, 1, 2, 3, 4, 5};

    EXPECT_SAME_TYPE(decltype(take3 | take1),
                     (seqan3::detail::combined_adaptor<take_adaptor_closure, take_adaptor_closure>));

    EXPECT_SAME_TYPE(decltype(take1(take3(vec))), decltype(std::views::take(std::views::take(vec, 3), 1)));
    EXPECT_SAME_TYPE(decltype(vec | (take3 | take1)), decltype(std::views::take(std::views::take(vec, 3), 1)));
    EXPECT_SAME_TYPE(decltype((take3 | take1)(vec)), decltype(std::views::take(std::views::take(vec, 3), 1)));

    EXPECT_RANGE_EQ(take1(take3(vec)), (std::vector<int>{0}));
    EXPECT_RANGE_EQ(vec | (take3 | take1), (std::vector<int>{0}));
    EXPECT_RANGE_EQ((take3 | take1)(vec), (std::vector<int>{0}));

    seqan3::detail::combined_adaptor adaptor = take3 | take1;
    EXPECT_SAME_TYPE(decltype(adaptor), (seqan3::detail::combined_adaptor<take_adaptor_closure, take_adaptor_closure>));

    EXPECT_SAME_TYPE(decltype(adaptor(vec)), decltype(std::views::take(std::views::take(vec, 3), 1)));
    EXPECT_SAME_TYPE(decltype(vec | adaptor), decltype(std::views::take(std::views::take(vec, 3), 1)));

    EXPECT_RANGE_EQ(adaptor(vec), (std::vector<int>{0}));
    EXPECT_RANGE_EQ(vec | adaptor, (std::vector<int>{0}));
}

TEST(adaptor_base, pipe_different_adaptor)
{
    constexpr take_adaptor_closure take1{1};
    constexpr drop_adaptor_closure drop3{3};
    std::vector<int> vec{0, 1, 2, 3, 4, 5};

    EXPECT_SAME_TYPE(decltype(drop3 | take1),
                     (seqan3::detail::combined_adaptor<drop_adaptor_closure, take_adaptor_closure>));

    EXPECT_SAME_TYPE(decltype(take1(drop3(vec))), decltype(std::views::take(std::views::drop(vec, 3), 1)));
    EXPECT_SAME_TYPE(decltype(vec | (drop3 | take1)), decltype(std::views::take(std::views::drop(vec, 3), 1)));
    EXPECT_SAME_TYPE(decltype((drop3 | take1)(vec)), decltype(std::views::take(std::views::drop(vec, 3), 1)));

    EXPECT_RANGE_EQ((drop3 | take1)(vec), (std::vector<int>{3}));
    EXPECT_RANGE_EQ(vec | (drop3 | take1), (std::vector<int>{3}));
    EXPECT_RANGE_EQ((drop3 | take1)(vec), (std::vector<int>{3}));

    seqan3::detail::combined_adaptor adaptor = drop3 | take1;
    EXPECT_SAME_TYPE(decltype(adaptor), (seqan3::detail::combined_adaptor<drop_adaptor_closure, take_adaptor_closure>));

    EXPECT_SAME_TYPE(decltype(adaptor(vec)), decltype(std::views::take(std::views::drop(vec, 3), 1)));
    EXPECT_SAME_TYPE(decltype(vec | adaptor), decltype(std::views::take(std::views::drop(vec, 3), 1)));

    EXPECT_RANGE_EQ(adaptor(vec), (std::vector<int>{3}));
    EXPECT_RANGE_EQ(vec | adaptor, (std::vector<int>{3}));
}

TEST(adaptor_base, pipe_left_non_seqan_adaptor)
{
    constexpr take_adaptor_closure take1{1};
    std::vector<int> vec{0, 1, 2, 3, 4, 5};

    EXPECT_SAME_TYPE(decltype(std::views::take(3) | take1),
                     (seqan3::detail::combined_adaptor<decltype(std::views::take(3)), take_adaptor_closure>));

    EXPECT_SAME_TYPE(decltype(std::views::take(3)(take1(vec))),
                     decltype(std::views::take(std::views::take(vec, 3), 1)));
    EXPECT_SAME_TYPE(decltype(vec | (std::views::take(3) | take1)),
                     decltype(std::views::take(std::views::take(vec, 3), 1)));
    EXPECT_SAME_TYPE(decltype((std::views::take(3) | take1)(vec)),
                     decltype(std::views::take(std::views::take(vec, 3), 1)));

    EXPECT_RANGE_EQ(std::views::take(3)(take1(vec)), (std::vector<int>{0}));
    EXPECT_RANGE_EQ(vec | (std::views::take(3) | take1), (std::vector<int>{0}));
    EXPECT_RANGE_EQ((std::views::take(3) | take1)(vec), (std::vector<int>{0}));

    seqan3::detail::combined_adaptor adaptor = std::views::take(3) | take1;
    EXPECT_SAME_TYPE(decltype(adaptor),
                     (seqan3::detail::combined_adaptor<decltype(std::views::take(3)), take_adaptor_closure>));

    EXPECT_SAME_TYPE(decltype(adaptor(vec)), decltype(std::views::take(std::views::take(vec, 3), 1)));
    EXPECT_SAME_TYPE(decltype(vec | adaptor), decltype(std::views::take(std::views::take(vec, 3), 1)));

    EXPECT_RANGE_EQ(adaptor(vec), (std::vector<int>{0}));
    EXPECT_RANGE_EQ(vec | adaptor, (std::vector<int>{0}));
}

TEST(adaptor_base, pipe_right_non_seqan_adaptor)
{
    constexpr take_adaptor_closure take1{1};
    std::vector<int> vec{0, 1, 2, 3, 4, 5};

    EXPECT_SAME_TYPE(decltype(take1 | std::views::take(3)),
                     (seqan3::detail::combined_adaptor<take_adaptor_closure, decltype(std::views::take(3))>));

    EXPECT_SAME_TYPE(decltype(take1(std::views::take(vec, 3))),
                     decltype(std::views::take(std::views::take(vec, 3), 1)));
    EXPECT_SAME_TYPE(decltype(vec | (take1 | std::views::take(3))),
                     decltype(std::views::take(std::views::take(vec, 3), 1)));
    EXPECT_SAME_TYPE(decltype((take1 | std::views::take(3))(vec)),
                     decltype(std::views::take(std::views::take(vec, 3), 1)));

    EXPECT_RANGE_EQ(take1(std::views::take(vec, 3)), (std::vector<int>{0}));
    EXPECT_RANGE_EQ(vec | (take1 | std::views::take(3)), (std::vector<int>{0}));
    EXPECT_RANGE_EQ((take1 | std::views::take(3))(vec), (std::vector<int>{0}));

    seqan3::detail::combined_adaptor adaptor = take1 | std::views::take(3);
    EXPECT_SAME_TYPE(decltype(adaptor),
                     (seqan3::detail::combined_adaptor<take_adaptor_closure, decltype(std::views::take(3))>));

    EXPECT_SAME_TYPE(decltype(adaptor(vec)), decltype(std::views::take(std::views::take(vec, 3), 1)));
    EXPECT_SAME_TYPE(decltype(vec | adaptor), decltype(std::views::take(std::views::take(vec, 3), 1)));

    EXPECT_RANGE_EQ(adaptor(vec), (std::vector<int>{0}));
    EXPECT_RANGE_EQ(vec | adaptor, (std::vector<int>{0}));
}

TEST(adaptor_base, rvalue_pipes)
{
    std::vector<int> vec{0, 1, 2, 3, 4, 5};

    auto take1_lvalue = take_adaptor_closure{1};
    auto take3_lvalue = std::views::take(3);

    auto take1_rvalue_take3_rvalue = take_adaptor_closure{1} | std::views::take(3);
    auto take1_rvalue_take3_lvalue = take_adaptor_closure{1} | take3_lvalue;
    auto take3_rvalue_take1_rvalue = std::views::take(3) | take_adaptor_closure{1};
    auto take3_rvalue_take1_lvalue = std::views::take(3) | take1_lvalue;

    EXPECT_SAME_TYPE(decltype(take1_rvalue_take3_rvalue),
                     (seqan3::detail::combined_adaptor<take_adaptor_closure, decltype(std::views::take(3))>));
    EXPECT_SAME_TYPE(decltype(take1_rvalue_take3_lvalue),
                     (seqan3::detail::combined_adaptor<take_adaptor_closure, decltype(std::views::take(3))>));
    EXPECT_SAME_TYPE(decltype(take3_rvalue_take1_rvalue),
                     (seqan3::detail::combined_adaptor<decltype(std::views::take(3)), take_adaptor_closure>));
    EXPECT_SAME_TYPE(decltype(take3_rvalue_take1_lvalue),
                     (seqan3::detail::combined_adaptor<decltype(std::views::take(3)), take_adaptor_closure>));

    EXPECT_RANGE_EQ(take1_rvalue_take3_rvalue(vec), (std::vector<int>{0}));
    EXPECT_RANGE_EQ(take1_rvalue_take3_lvalue(vec), (std::vector<int>{0}));
    EXPECT_RANGE_EQ(take3_rvalue_take1_rvalue(vec), (std::vector<int>{0}));
    EXPECT_RANGE_EQ(take3_rvalue_take1_lvalue(vec), (std::vector<int>{0}));
}
