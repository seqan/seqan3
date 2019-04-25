// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/range/container/constexpr_string.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/std/ranges>

#include <utility>
#include <string>

using namespace std::literals;
using namespace seqan3;

// standard construction.
TEST(constexpr_string, standard_construction)
{
    EXPECT_TRUE((std::is_default_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_trivially_default_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_nothrow_default_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_copy_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_trivially_copy_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_nothrow_copy_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_move_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_trivially_move_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_nothrow_move_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_copy_assignable_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_trivially_copy_assignable_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_nothrow_copy_assignable_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_move_assignable_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_trivially_move_assignable_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_nothrow_move_assignable_v<constexpr_string<4>>));
}

TEST(constexpr_string, Container)
{
    EXPECT_TRUE(Container<constexpr_string<4>>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<constexpr_string<4>>);
}

// construction from literal.
TEST(constexpr_string, construct_from_literal)
{
     EXPECT_TRUE((std::is_same_v<constexpr_string<5>, decltype(constexpr_string{"hello"})>));
}

// construction from char.
TEST(constexpr_string, construct_from_char)
{
     EXPECT_TRUE((std::is_same_v<constexpr_string<1>, decltype(constexpr_string{'h'})>));
}

// construction from array.
TEST(constexpr_string, construct_from_array)
{
     EXPECT_TRUE((std::is_same_v<constexpr_string<5>,
                                 decltype(constexpr_string{std::array{'h','e','l','l','o'}})>));
}

TEST(constexpr_string, size)
{
    constexpr_string em{"hello"};

    EXPECT_EQ(em.size(), 5u);
    constexpr auto size = em.size();
    EXPECT_EQ(size, 5u);
}

TEST(constexpr_string, max_size)
{
    constexpr_string em{"hello"};

    EXPECT_EQ(em.max_size(), 5u);
    constexpr auto size = em.size();
    EXPECT_EQ(size, 5u);
}

TEST(constexpr_string, c_str)
{
    {
        constexpr_string em{"hello"};
        EXPECT_EQ(std::string{em.c_str()}, "hello"s);
    }

    {
        constexpr_string em{'x'};
        EXPECT_EQ(std::string{em.c_str()}, "x"s);
    }
}

TEST(constexpr_string, string)
{
    constexpr_string em{"hello"};
    EXPECT_EQ(em.string(), "hello"s);  // explicit
}

TEST(constexpr_string, concat)
{
    {
        constexpr_string em = constexpr_string{"hello"} +
                              constexpr_string{' '} +
                              constexpr_string{"world"};
        constexpr auto size = em.size();

        EXPECT_EQ(size, 11u);
        EXPECT_EQ(em.string(), "hello world"s);
    }

    {
        static constexpr char const a[] = "hello";
        static constexpr char const b[] = " ";
        static constexpr char const c[] = "world";
        auto em = constexpr_string{a} + constexpr_string{b} + constexpr_string{c};
        constexpr auto size = em.size();

        EXPECT_EQ(size, 11u);
        EXPECT_EQ(em.string(), "hello world"s);
    }
}

TEST(constexpr_string, begin)
{
    constexpr_string s{"hello"};
    EXPECT_EQ(*s.begin(), 'h');
    EXPECT_TRUE((std::is_same_v<decltype(s)::iterator, decltype(s.begin())>));

    constexpr_string<5> const cs{s};
    EXPECT_EQ(*cs.begin(), 'h');
    EXPECT_TRUE((std::is_same_v<decltype(cs)::const_iterator, decltype(cs.begin())>));
}

TEST(constexpr_string, cbegin)
{
    constexpr_string s{"hello"};
    EXPECT_EQ(*s.begin(), 'h');
    EXPECT_TRUE((std::is_same_v<decltype(s)::const_iterator, decltype(s.cbegin())>));
}

TEST(constexpr_string, end)
{
    constexpr_string s{"hello"};
    EXPECT_EQ(*(s.end() - 1), 'o');
    EXPECT_TRUE((std::is_same_v<decltype(s)::iterator, decltype(s.end())>));

    constexpr_string<5> const cs{s};
    static_assert(std::is_same_v<constexpr_string<5> const, decltype(cs)>);
    EXPECT_EQ(*(cs.end() - 1), 'o');
    EXPECT_TRUE((std::is_same_v<decltype(cs)::const_iterator, decltype(cs.end())>));
}

TEST(constexpr_string, cend)
{
    constexpr_string s{"hello"};
    EXPECT_EQ(*(s.cend() - 1), 'o');
    EXPECT_TRUE((std::is_same_v<decltype(s)::const_iterator, decltype(s.cend())>));
}

TEST(constexpr_string, swap)
{
    constexpr_string s1{"hello"};
    constexpr_string s2{"olleh"};
    {  // global function.
        std::swap(s1, s2);
        EXPECT_EQ(s1, constexpr_string{"olleh"});
        EXPECT_EQ(s2, constexpr_string{"hello"});
    }

    { // member function
        s1.swap(s2);
        EXPECT_EQ(s1, constexpr_string{"hello"});
        EXPECT_EQ(s2, constexpr_string{"olleh"});
    }
}

TEST(constexpr_string, equality)
{
    constexpr bool cmp1 = constexpr_string{"hello"} == constexpr_string{"hello"};
    constexpr bool cmp2 = constexpr_string{"hello"} == constexpr_string{"hell"};
    constexpr bool cmp3 = constexpr_string{"hell"}  == constexpr_string{"hello"};
    constexpr bool cmp4 = constexpr_string{"hella"} == constexpr_string{"hello"};

    EXPECT_TRUE(cmp1);
    EXPECT_FALSE(cmp2);
    EXPECT_FALSE(cmp3);
    EXPECT_FALSE(cmp4);
}

TEST(constexpr_string, inequality)
{
    constexpr bool cmp1 = constexpr_string{"hello"} != constexpr_string{"hello"};
    constexpr bool cmp2 = constexpr_string{"hello"} != constexpr_string{"hell"};
    constexpr bool cmp3 = constexpr_string{"hell"}  != constexpr_string{"hello"};
    constexpr bool cmp4 = constexpr_string{"hella"} != constexpr_string{"hello"};

    EXPECT_FALSE(cmp1);
    EXPECT_TRUE(cmp2);
    EXPECT_TRUE(cmp3);
    EXPECT_TRUE(cmp4);
}

TEST(constexpr_string, less)
{
    constexpr bool cmp1 = constexpr_string{"hello"} < constexpr_string{"hello"};
    constexpr bool cmp2 = constexpr_string{"hello"} < constexpr_string{"hell"};
    constexpr bool cmp3 = constexpr_string{"hell"}  < constexpr_string{"hello"};
    constexpr bool cmp4 = constexpr_string{"hella"} < constexpr_string{"hello"};

    EXPECT_FALSE(cmp1);
    EXPECT_FALSE(cmp2);
    EXPECT_TRUE(cmp3);
    EXPECT_TRUE(cmp4);
}

TEST(constexpr_string, less_equal)
{
    constexpr bool cmp1 = constexpr_string{"hello"} <= constexpr_string{"hello"};
    constexpr bool cmp2 = constexpr_string{"hello"} <= constexpr_string{"hell"};
    constexpr bool cmp3 = constexpr_string{"hell"}  <= constexpr_string{"hello"};
    constexpr bool cmp4 = constexpr_string{"hella"} <= constexpr_string{"hello"};

    EXPECT_TRUE(cmp1);
    EXPECT_FALSE(cmp2);
    EXPECT_TRUE(cmp3);
    EXPECT_TRUE(cmp4);
}

TEST(constexpr_string, greater)
{
    constexpr bool cmp1 = constexpr_string{"hello"} > constexpr_string{"hello"};
    constexpr bool cmp2 = constexpr_string{"hello"} > constexpr_string{"hell"};
    constexpr bool cmp3 = constexpr_string{"hell"}  > constexpr_string{"hello"};
    constexpr bool cmp4 = constexpr_string{"hella"} > constexpr_string{"hello"};

    EXPECT_FALSE(cmp1);
    EXPECT_TRUE(cmp2);
    EXPECT_FALSE(cmp3);
    EXPECT_FALSE(cmp4);
}

TEST(constexpr_string, greater_equal)
{
    constexpr bool cmp1 = constexpr_string{"hello"} >= constexpr_string{"hello"};
    constexpr bool cmp2 = constexpr_string{"hello"} >= constexpr_string{"hell"};
    constexpr bool cmp3 = constexpr_string{"hell"}  >= constexpr_string{"hello"};
    constexpr bool cmp4 = constexpr_string{"hella"} >= constexpr_string{"hello"};

    EXPECT_TRUE(cmp1);
    EXPECT_TRUE(cmp2);
    EXPECT_FALSE(cmp3);
    EXPECT_FALSE(cmp4);
}

template <std::size_t N>
constexpr constexpr_string<N> fill_constexpr_string(constexpr_string<N> s, char const val)
{
    for (auto && c : s)
    {
        c = val;
    }
    return s;
}

TEST(constexpr_string, compile_time_fill)
{
    constexpr bool cmp = fill_constexpr_string(constexpr_string<4>{}, 'x') == constexpr_string{"xxxx"};
    EXPECT_TRUE(cmp);
}
